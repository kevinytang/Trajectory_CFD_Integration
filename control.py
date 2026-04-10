"""
Coupled Trajectory CFD Controller
==========================================================
Author: Kevin Tang

Usage (single manual run — original behaviour):
  cd Trajectory_CFD_Integration/
  python control.py

Usage (pipeline / multi-case run):
  python control.py --workdir runs/case_001 --config runs/case_001/config.nml \\
                    --mesh-dir fun3d_meshes/cone_59deg

Instructions:
1. Check if all the paths are correct.
2. Adjust input parameters as needed.
3. Make sure T_step matches t_end in the fortran code.
4. Enter the ESP shell environment.
5. Use the virtual environment (.venv) in the current directory.
"""

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

########################################
##         Helper Functions           ##
########################################

def readNML(CFG):
    params = f90nml.read(CFG)
    traj = params["current_states"]
    cfd  = params["cfd_variables"]
    return {
        "step":        traj.get("step"),
        "time":        traj.get("time"),
        "v":           traj.get("v"),
        "gamma":       traj.get("gamma"),
        "psi":         traj.get("psi"),
        "alt":         traj.get("alt"),
        "lon":         traj.get("lon"),
        "lat":         traj.get("lat"),
        "cl":          traj.get("cl"),
        "cd":          traj.get("cd"),
        "mach":        cfd.get("mach"),
        "temperature": cfd.get("temperature"),
        "pressure":    cfd.get("pressure"),
        "density":     cfd.get("density"),
        "alpha":       cfd.get("alpha"),
    }


def writeNML(CFG, step, time, v, gamma, psi, alt, lon, lat,
             cl, cd, mach, temperature, pressure, density, alpha):
    params = f90nml.read(CFG)
    traj   = params["current_states"]
    cfd    = params["cfd_variables"]

    traj["step"] = step;  traj["time"]  = time
    traj["v"]    = v;     traj["gamma"] = gamma;  traj["psi"] = psi
    traj["alt"]  = alt;   traj["lon"]   = lon;    traj["lat"] = lat
    traj["cl"]   = cl;    traj["cd"]    = cd

    cfd["mach"]        = mach
    cfd["temperature"] = temperature
    cfd["pressure"]    = pressure
    cfd["density"]     = density
    cfd["alpha"]       = alpha

    with open(CFG, "w") as f:
        f90nml.write(params, f)
    return params


def writeCSV(HIST, step, time, v, gamma, psi, alt, lon, lat,
             cl, cd, mach, temperature, pressure, density, alpha):
    fieldnames = ["step", "time", "v", "gamma", "psi", "alt", "lon", "lat",
                  "cl", "cd", "mach", "temperature", "pressure", "density", "alpha"]
    data = {
        "step": step, "time": time, "v": v, "gamma": gamma, "psi": psi,
        "alt": alt, "lon": lon, "lat": lat, "cl": cl, "cd": cd,
        "mach": mach, "temperature": temperature, "pressure": pressure,
        "density": density, "alpha": alpha,
    }
    write_header = not HIST.exists()
    with open(HIST, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(data)
    return data


########################################
##         CLI argument parsing       ##
########################################

def _parse_args():
    parser = argparse.ArgumentParser(description="Coupled Trajectory-CFD Controller")
    parser.add_argument("--workdir",  default=None,
                        help="Working directory for this run (default: script's directory). "
                             "All relative paths are resolved from here.")
    parser.add_argument("--config",   default=None,
                        help="Path to config.nml (default: <workdir>/config.nml)")
    parser.add_argument("--mesh-dir", default=None,
                        help="Path to pre-built FUN3D mesh directory. "
                             "Passed through to run_fun3d.py to skip mesh generation.")
    return parser.parse_args()


########################################
##         Main execution             ##
########################################

def main():
    args = _parse_args()

    # ── Resolve working directory ────────────────────────────────────────────
    script_dir = Path(__file__).parent.resolve()
    workdir    = Path(args.workdir).resolve() if args.workdir else script_dir
    workdir.mkdir(parents=True, exist_ok=True)

    # ── Resolve file paths relative to workdir ───────────────────────────────
    CFG    = Path(args.config).resolve() if args.config else workdir / "config.nml"
    HIST   = workdir / "history.csv"

    # Trajectory executable and CFD script are always relative to the submodule dir
    TRAJ   = script_dir / "Reentry_3DOF_NonPlanar" / "reentry.exe"
    RUNCFD = script_dir / "fun3D_Solver" / "run_fun3d.py"

    # Mesh directory for FUN3D (optional — enables mesh reuse)
    mesh_dir = Path(args.mesh_dir).resolve() if args.mesh_dir else None

    # reentry.exe opens '../config.nml' relative to its CWD (hardcoded in Fortran).
    # Run it from workdir/reentry_work/ so that '../config.nml' resolves to
    # workdir/config.nml (the case-specific config) rather than the submodule template.
    reentry_workdir = workdir / "reentry_work"
    reentry_workdir.mkdir(parents=True, exist_ok=True)

    print(f"[control] workdir  : {workdir}")
    print(f"[control] config   : {CFG}")
    print(f"[control] mesh-dir : {mesh_dir or 'None (mesh will be generated)'}")

    if not CFG.exists():
        print(f"[control] ERROR: config.nml not found at {CFG}", file=sys.stderr)
        sys.exit(1)

    # ── Read all settings from config.nml ───────────────────────────────────
    params = f90nml.read(CFG)

    init  = params["initial_conditions"]
    v     = init["v0"];    gamma = init["gamma0"]; psi   = init["psi0"]
    alt   = init["alt0"];  lon   = init["lon0"];   lat   = init["lat0"]
    alpha = init["alpha0"]

    ctrl   = params["control_settings"]
    t_step = ctrl["t_step"]
    t_end  = ctrl["t_end_ctrl"]
    tol    = ctrl["tol"]

    step = 0;  time = 0
    cl = cd = mach = temperature = pressure = density = 0.0

    # Build CFD subprocess command
    cfd_cmd = ["python", str(RUNCFD), "--config", str(CFG)]
    if mesh_dir:
        cfd_cmd += ["--mesh-dir", str(mesh_dir)]

    # ── Write initial namelist state ─────────────────────────────────────────
    writeNML(CFG, step, time, v, gamma, psi, alt, lon, lat,
             cl, cd, mach, temperature, pressure, density, alpha)

    # ── Initial trajectory step → compute Mach / Temperature ─────────────────
    subprocess.run([str(TRAJ)], cwd=str(reentry_workdir), check=True)
    state = readNML(CFG)

    # ── Initial CFD run → get first Cl/Cd ───────────────────────────────────
    subprocess.run(cfd_cmd, check=True)
    state      = readNML(CFG)
    CFD_input  = np.array([state["mach"], state["temperature"], state["alpha"]])

    writeCSV(HIST, state["step"], state["time"], state["v"], state["gamma"],
             state["psi"], state["alt"], state["lon"], state["lat"],
             state["cl"], state["cd"], state["mach"], state["temperature"],
             state["pressure"], state["density"], state["alpha"])

    # ── Main propagation loop ─────────────────────────────────────────────────
    while time < t_end:

        subprocess.run([str(TRAJ)], cwd=str(reentry_workdir), check=True)
        state = readNML(CFG)

        if state["alt"] <= 0:
            print(f"[LANDED] t={state['time']:.1f}s | alt={state['alt']:.1f}m")
            writeCSV(HIST, state["step"], state["time"], state["v"], state["gamma"],
                     state["psi"], state["alt"], state["lon"], state["lat"],
                     state["cl"], state["cd"], state["mach"], state["temperature"],
                     state["pressure"], state["density"], state["alpha"])
            break

        CFD_input_new = np.array([state["mach"], state["temperature"], state["alpha"]])
        denom         = np.where(CFD_input != 0, CFD_input, 1)
        per_diff      = np.abs(CFD_input_new - CFD_input) / denom * 100
        outside_tol   = per_diff > tol

        if np.any(outside_tol):
            subprocess.run(cfd_cmd, check=True)
            CFD_input = CFD_input_new

        writeCSV(HIST, state["step"], state["time"], state["v"], state["gamma"],
                 state["psi"], state["alt"], state["lon"], state["lat"],
                 state["cl"], state["cd"], state["mach"], state["temperature"],
                 state["pressure"], state["density"], state["alpha"])

        step += 1
        time += t_step

        print(f"[Step {int(state['step']):04d}] t={state['time']:.1f}s | "
              f"max ΔCFD={np.max(per_diff):.2f}% | CFD re-run={np.any(outside_tol)}")

    print(f"[control] Done. History written to {HIST}")


if __name__ == "__main__":
    main()
