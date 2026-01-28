
"""
Coupled Trajectory CFD Controller
==========================================================
Author: Kevin Tang
Instructions:
1. Check if all the paths are correct.
2. Adjust input parameters as needed.
3. Make sure T_step matches t_end in the fortran code.
4. Enter the ESP shell environment.
5. Use the virtual environment (.venv) in the current directory.
"""

import csv
import subprocess
import f90nml
import math
import numpy as np
from pathlib import Path

########################################
##         Helper Functions           ##
########################################

# Read from the namelist
def readNML(CFG):

    # Get files
    params = f90nml.read(CFG)

    # Add access to each sections
    traj = params["current_states"]
    cfd = params["cfd_variables"]

    # Update each variables
    step = traj["step"]
    time = traj["time"]
    v = traj["v"]
    gamma = traj["gamma"]
    psi = traj["psi"]
    alt = traj["alt"]
    lon = traj["lon"]
    lat = traj["lat"]
    cl = traj["cl"]
    cd = traj["cd"]

    mach = cfd["mach"]
    temperature = cfd["temperature"]
    alpha = cfd["alpha"]

    return {
        "step": step, "time": time, "v": v, "gamma": gamma, "psi": psi, "alt": alt, "lon": lon, "lat": lat,
        "cl": cl, "cd": cd, "mach": mach, "temperature": temperature, "alpha": alpha
    }

# Write to the namelist
def writeNML(CFG, step, time, v, gamma, psi, alt, lon, lat, cl, cd, mach, temperature, alpha):

    # Get files
    params = f90nml.read(CFG)

    # Add access to each sections
    traj = params["current_states"]
    cfd = params["cfd_variables"]

    # Update each variables
    traj["step"] = step
    traj["time"] = time
    traj["v"] = v
    traj["gamma"] = gamma
    traj["psi"] = psi
    traj["alt"] = alt
    traj["lon"] = lon
    traj["lat"] = lat
    traj["cl"] = cl
    traj["cd"] = cd

    cfd["mach"] = mach
    cfd["temperature"] = temperature
    cfd["alpha"] = alpha

    # Write the initial condition
    with open(CFG, "w") as f:
        f90nml.write(params, f)

    return params

# Write to csv file
def writeCSV(HIST, step, time, v, gamma, psi, alt, lon, lat, cl, cd, mach, temperature, alpha):

    # Define fieldnames once
    fieldnames = ["step", "time", "v", "gamma", "psi", "alt", "lon", "lat", "cl", "cd", "mach", "temperature",
                  "alpha"]

    # Data in one row
    data = {
        "step": step,
        "time": time,
        "v": v,
        "gamma": gamma,
        "psi": psi,
        "alt": alt,
        "lon": lon,
        "lat": lat,
        "cl": cl,
        "cd": cd,
        "mach": mach,
        "temperature": temperature,
        "alpha": alpha
    }

    # Write the data
    write_header = not HIST.exists()
    with open(HIST, "a", newline = "") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(data)

    return data

########################################
##         Initial Variables          ##
########################################
# Paths (adjust to your folder layout)
CFG   = Path("config.nml")
RUNCFD= Path("fun3D_Solver/run_075_70deg.py")
TRAJ  = Path("Reentry_3DOF_NonPlanar/reentry.exe")
HIST = Path("history.csv")

# Read all settings from config.nml
params = f90nml.read(CFG)

# Initial conditions
init = params["initial_conditions"]
v = init["v0"]
gamma = init["gamma0"]
psi = init["psi0"]
alt = init["alt0"]
lon = init["lon0"]
lat = init["lat0"]

# Control settings
ctrl = params["control_settings"]
t_step = ctrl["t_step"]
t_end = ctrl["t_end_ctrl"]
tol = ctrl["tol"] # Percent difference threshold to call CFD

# Initial condition
step = 0
time = 0

# Placeholders (will be computed by Fortran after first run)
cl = 0.0
cd = 0.0
mach = 0.0
temperature = 0.0
alpha = gamma  # Initial angle of attack = flight path angle

########################################
##       Write the Initial Input      ##
########################################
writeNML(CFG, step, time, v, gamma, psi, alt, lon, lat, cl, cd, mach, temperature, alpha)

########################################
##         Running Propagation        ##
########################################

# Run trajectory propagator ONCE with one timestep to compute initial Mach & Temperature
subprocess.run(['./' + TRAJ.name], cwd=TRAJ.parent)

# Read back the computed Mach and Temperature
state = readNML(CFG)

# Run CFD ONCE to get the coefficients
subprocess.run(["python", RUNCFD.name], cwd=RUNCFD.parent)

# Record CFD input and write the initial state to the csv file
state = readNML(CFG)
CFD_input = np.array([state["mach"], state["temperature"], state["alpha"]])
writeCSV(HIST, state["step"], state["time"], state["v"], state["gamma"], state["psi"], state["alt"], state["lon"],
         state["lat"], state["cl"], state["cd"], state["mach"], state["temperature"], state["alpha"])

# While loop
while time < t_end:

    # Run the trajectory propagator for n seconds
    subprocess.run(['./' + TRAJ.name], cwd=TRAJ.parent)

    # Get the current state
    state = readNML(CFG)

    # Read the CFD inputs from the namelist file and find the percent difference
    CFD_input_new = np.array([state["mach"], state["temperature"], state["alpha"]])
    denom = np.where(CFD_input != 0, CFD_input, 1) # guard against 0
    per_diff = np.abs(CFD_input_new - CFD_input) / denom * 100

    # Boolean flag: check if any percent difference exceed the tol
    outside_tol = per_diff > tol

    # If: CFD inputs are too different from the last recorded input
    if np.any(outside_tol):

        # Run CFD again to update Cl and Cd
        subprocess.run(["python", RUNCFD.name], cwd=RUNCFD.parent)

        # Save the CFD input for future reference
        CFD_input = CFD_input_new

    # Write the current state to the csv file
    writeCSV(HIST, state["step"], state["time"], state["v"], state["gamma"], state["psi"], state["alt"], state["lon"],
             state["lat"], state["cl"], state["cd"], state["mach"], state["temperature"], state["alpha"])

    # Update and prepare for the next loop
    step = step + 1
    time = time + t_step

    # Output CFD re-run
    print(f"[Step {int(state['step']):04d}] t={state['time']:.1f}s | max ΔCFD={np.max(per_diff):.2f}% | CFD re-run={np.any(outside_tol)}")
