"""
run_fun3d.py
============
FUN3D CFD driver via pyCAPS (ESP Engineering Sketch Pad).

Reads trajectory state from config.nml, configures FUN3D for inviscid
hypersonic flow, runs the solver, and writes updated CL/CD back to config.nml.

Mesh reuse:
  If --mesh-dir is provided and contains a mesh_ready.flag, the AFLR4/AFLR3
  mesh generation steps are skipped and the pre-built mesh is loaded directly.
  This avoids regenerating the FUN3D 3D volume mesh on every CFD call.

  NOTE: The FUN3D 3D volume mesh is completely separate from the FEAR 2D Gmsh
  mesh.  Do not confuse the two.

Usage:
  # Standard (generates mesh each call):
  python run_fun3d.py --config ../config.nml

  # With pre-built mesh (pipeline mode):
  python run_fun3d.py --config runs/case_001/config.nml \\
                      --mesh-dir fun3d_meshes/cone_59deg
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

import f90nml
import numpy as np
import pyCAPS

FLAG_FILE = "mesh_ready.flag"


def fbool(x: bool) -> str:
    return ".true." if bool(x) else ".false."


def _parse_args():
    parser = argparse.ArgumentParser(description="FUN3D CFD driver (pyCAPS)")
    parser.add_argument("--config",   default=None,
                        help="Path to config.nml "
                             "(default: ../config.nml relative to this script)")
    parser.add_argument("--mesh-dir", default=None,
                        help="Pre-built FUN3D mesh directory. "
                             "If valid and mesh_ready.flag exists, skip AFLR4/AFLR3.")
    return parser.parse_args()


def main():
    args    = _parse_args()
    cfd_dir = Path(__file__).parent.resolve()

    # ── Resolve config.nml path ──────────────────────────────────────────────
    if args.config:
        config_path = Path(args.config).resolve()
    else:
        config_path = (cfd_dir / ".." / "config.nml").resolve()

    if not config_path.exists():
        print(f"[run_fun3d] ERROR: config.nml not found at {config_path}", file=sys.stderr)
        sys.exit(1)

    # ── Load parameters ──────────────────────────────────────────────────────
    params = f90nml.read(str(config_path))
    cfdvar = f90nml.read(str(cfd_dir / "cfd_control.nml"))

    traj = params["current_states"]
    cfd  = params["cfd_variables"]
    var  = cfdvar["cfd_var"]

    # ── Determine mesh mode ──────────────────────────────────────────────────
    mesh_dir   = Path(args.mesh_dir).resolve() if args.mesh_dir else None
    reuse_mesh = (mesh_dir is not None and
                  (mesh_dir / FLAG_FILE).exists())

    if reuse_mesh:
        print(f"[run_fun3d] Reusing pre-built FUN3D mesh from {mesh_dir}")
    else:
        print("[run_fun3d] Generating FUN3D mesh (AFLR4 + AFLR3)")

    # ── pyCAPS problem name: unique per case to avoid collisions ─────────────
    # Derive a unique suffix from the config path so parallel runs don't clash
    case_suffix = config_path.parent.name.replace("/", "_")
    prob_name   = f"reentryCFD_{case_suffix}"

    filename = var["filename"]
    csm_path = cfd_dir / filename
    if not csm_path.exists():
        print(f"[run_fun3d] ERROR: CSM file not found: {csm_path}", file=sys.stderr)
        sys.exit(1)

    caps_prob = pyCAPS.Problem(problemName=prob_name,
                               capsFile=str(csm_path),
                               outLevel=1)

    if reuse_mesh:
        # ── Load pre-built mesh directly ─────────────────────────────────────
        # Link a dummy AFLR3 AIM that reads the existing .b8.ugrid
        aflr3_name = "aflr3"
        aflr3 = caps_prob.analysis.create(aim="aflr3AIM", name=aflr3_name)

        ugrid_files = list(mesh_dir.glob("*.b8.ugrid")) + list(mesh_dir.glob("*.lb8.ugrid"))
        if not ugrid_files:
            print(f"[run_fun3d] WARNING: No .b8.ugrid in {mesh_dir}; falling back to remesh")
            reuse_mesh = False
        else:
            # Copy pre-built mesh into pyCAPS scratch dir so FUN3D AIM can find it
            scratch = Path(caps_prob.analysisDir(aflr3_name))
            scratch.mkdir(parents=True, exist_ok=True)
            for uf in ugrid_files:
                shutil.copy2(uf, scratch / uf.name)
            # Also copy .mapbc if present
            for mf in mesh_dir.glob("*.mapbc"):
                shutil.copy2(mf, scratch / mf.name)
            aflr3.runAnalysis()

    if not reuse_mesh:
        # ── Build mesh from scratch ───────────────────────────────────────────
        aflr4 = caps_prob.analysis.create(aim="aflr4AIM", name="aflr4")
        aflr4.input.Mesh_Length_Factor = var.get("Mesh_Length_Factor")
        aflr4.input.max_scale          = var.get("max_scale")
        aflr4.input.ideal_min_scale    = var.get("min_scale")
        aflr4.input.ff_cdfr            = var.get("ff_cdfr")
        aflr4.input.Mesh_Sizing = {
            "blunt":    {"edgeWeight":  var.get("edgeWeight"),
                         "scaleFactor": var.get("blunt_scaleFactor")},
            "Farfield": {"bcType":      "Farfield",
                         "scaleFactor": var.get("farfield_scaleFactor")},
        }
        print("==> Running AFLR4 …")
        aflr4.runAnalysis()

        aflr3 = caps_prob.analysis.create(aim="aflr3AIM", name="aflr3")
        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])
        aflr3.input.Mesh_Sizing = {
            "blunt":    {"bcType": "Inviscid"},
            "Farfield": {"bcType": "Farfield"},
        }
        print("==> Running AFLR3 …")
        aflr3.runAnalysis()

    # ── FUN3D AIM ────────────────────────────────────────────────────────────
    print("==> Configuring FUN3D AIM …")
    fun3d = caps_prob.analysis.create(aim="fun3dAIM", name="fun3d")

    proj_name = var.get("Proj_Name", "fun3d_run")
    fun3d.input.Proj_Name     = proj_name
    fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])
    fun3d.input.Use_Python_NML = var.get("Use_Python_NML", True)
    fun3d.input.Overwrite_NML  = var.get("Overwrite_NML",  True)

    nproc = int(os.environ.get("FUN3D_MPI_PROCS", str(var.get("np", 10))))

    fun3d.input.Num_Iter          = int(var.get("Num_Iter"))
    fun3d.input.CFL_Schedule      = var.get("CFL_Schedule")
    fun3d.input.CFL_Schedule_Iter = var.get("CFL_Schedule_Iter")
    fun3d.input.Restart_Read      = var.get("Restart_Read")

    fun3d.input.Boundary_Condition = {
        "blunt":    {"bcType": "Inviscid"},
        "Farfield": {"bcType": "Freestream",
                     "machNumber": cfd.get("mach")},
    }

    fun3d.preAnalysis()

    # ── Patch fun3d.nml with extra blocks ───────────────────────────────────
    r      = 0.4   # farfield sphere radius [m]
    A_ref  = np.pi * r**2
    var["Characteristic_Surface"] = float(A_ref)

    nml_path = os.path.join(fun3d.analysisDir, "fun3d.nml")
    with open(nml_path, "a") as f:
        f.write("&reference_physical_properties\n")
        f.write(f"  dim_input_type = '{var.get('dim_input_type')}'\n")
        f.write(f"  temperature_units = '{var.get('temperature_units')}'\n")
        f.write(f"  gridlength_conversion = {float(var.get('gridlength_conversion'))}\n")
        f.write(f"  mach_number = {float(cfd.get('mach'))}\n")
        f.write(f"  temperature = {float(cfd.get('temperature'))}\n")
        f.write(f"  angle_of_attack = {float(cfd.get('alpha'))}\n")
        f.write(f"  angle_of_yaw = {float(var.get('Beta'))}\n")
        f.write("/\n\n")

        f.write("&governing_equations\n")
        f.write(f"  eqn_type = '{var.get('Equation_Type')}'\n")
        f.write(f"  viscous_terms = '{var.get('Viscous')}'\n")
        f.write("/\n\n")

        f.write("&inviscid_flux_method\n")
        f.write(f"  flux_construction = '{var['Flux_Construction']}'\n")
        f.write(f"  flux_limiter = '{var['Flux_Limiter']}'\n")
        f.write(f"  first_order_iterations = {int(var['First_Order_Iterations'])}\n")
        f.write(f"  freeze_limiter_iteration = {int(var['Freeze_Limiter_Iteration'])}\n")
        f.write("/\n\n")

        f.write("&nonlinear_solver_parameters\n")
        f.write(f"  schedule_iteration = "
                f"{var['CFL_Schedule_Iter'][0]} {var['CFL_Schedule_Iter'][1]}\n")
        f.write(f"  schedule_cfl = "
                f"{var['CFL_Schedule'][0]} {var['CFL_Schedule'][1]}\n")
        f.write("/\n\n")

        f.write("&force_moment_integ_properties\n")
        f.write(f"  area_reference = {var['Characteristic_Surface']}\n")
        f.write("/\n\n")

        f.write("&linear_solver_control\n")
        f.write("  linear_projection = .true.\n")
        f.write("/\n\n")

    # ── Run FUN3D ────────────────────────────────────────────────────────────
    nodet = var.get("nodet_Dir")
    cmd   = (f"mpirun -np {nproc} {nodet} "
             f"--animation_freq -1 --volume_animation_freq -1")
    freeze = var.get("Freeze_Limiter")
    if freeze:
        cmd += f" --freeze_limiter {freeze}"

    print(f"[run_fun3d] Running: {cmd}")
    fun3d.system(cmd)
    fun3d.postAnalysis()

    # ── Extract CL / CD and write back to config.nml ─────────────────────────
    force_file = os.path.join(fun3d.analysisDir, f"{proj_name}.forces")
    CL = CD = None

    if os.path.exists(force_file):
        with open(force_file) as f:
            for line in reversed(f.readlines()):
                if line.strip().startswith("Cl") and "Cd" in line:
                    parts = line.strip().replace("=", " ").split()
                    CL = float(parts[1])
                    CD = float(parts[3])
                    break

    if CL is not None and CD is not None:
        fresh = f90nml.read(str(config_path))
        fresh["current_states"]["cl"] = CL
        fresh["current_states"]["cd"] = CD
        with open(config_path, "w") as f:
            f90nml.write(fresh, f)
        print(f"[run_fun3d] CL={CL:.6f}  CD={CD:.6f} → {config_path}")
    else:
        print(f"[run_fun3d] WARNING: Could not extract CL/CD from {force_file}")


if __name__ == "__main__":
    main()
