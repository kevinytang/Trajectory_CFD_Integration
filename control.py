
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
    traj = params["trajectory_input"]
    cfd = params["cfd_input"]

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
    traj = params["trajectory_input"]
    cfd = params["cfd_input"]

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


def earth_atm_nasa(alt_m, velocity = None):
    """
    Earth Standard Atmosphere model (NASA)
    ------------------------------------------------------------
    Computes temperature [°C], pressure [kPa], density [kg/m³],
    and speed of sound [m/s] at a given geometric altitude [m].
    Optionally computes Mach number if velocity [m/s] is given.

    Reference:
        https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    """

    # Constants
    T_ref = 273.15  # Reference temperature [K]
    gamma = 1.4  # Ratio of specific heats
    R = 287.0  # Specific gas constant [J/kg/K]

    # Determine temperature and pressure by layer
    if alt_m > 25000.0:
        # Upper stratosphere
        T_C = -131.21 + 0.00299 * alt_m
        P_kPa = 2.488 * ((T_C + 273.1) / 216.6) ** (-11.388)

    elif 11000.0 <= alt_m <= 25000.0:
        # Lower stratosphere
        T_C = -56.46
        P_kPa = 22.65 * math.exp(1.73 - 0.000157 * alt_m)

    else:
        # Troposphere
        T_C = 15.04 - 0.00649 * alt_m
        P_kPa = 101.29 * ((T_C + 273.1) / 288.08) ** 5.256

    # Speed of sound (CPG assumption)
    T_K = T_C + T_ref
    a = math.sqrt(gamma * R * T_K)

    # Density
    rho = P_kPa / (0.2869 * (T_C + 273.1))

    result = {
        "T_K": T_K,
        "T_C": T_C
    }

    # Optional Mach number
    if velocity is not None:
        result["Mach"] = velocity / a

    return result

########################################
##         Initial Variables          ##
########################################
# Paths (adjust to your folder layout)
CFG   = Path("config.nml")
RUNCFD= Path("fun3D_Solver/run_075_70deg.py")
TRAJ  = Path("Reentry_3DOF_NonPlanar/reentry.exe")
HIST = Path("history.csv")

# Initial condition
step = 0
time = 0
v = 11032.05312 # velocity [m/s]
gamma = -6.48 # flight path angle [deg]
psi = 39.82 # heading angle [deg]
alt = 121920 # altitude [m]
lon = 171.96 # longitude [deg]
lat = -3.19 # latitude [deg]

# Initial derived quantity for CFD
atm = earth_atm_nasa(alt, velocity=v)
temperature = atm["T_K"]
mach = atm["Mach"]
alpha = gamma

# Placeholder for aerodynamic coefficients (will be updated after CFD)
cl = 0 # lift coefficient
cd = 0 # drag coefficient

# Propagator control
t_step = 5 # time of each step for the propagator (also change "t_end" in fortran code for consistency)
t_end = 540 # total simulation time
tol = 10 # threshold on calling CFD in percent

########################################
##       Write the Initial Input      ##
########################################
writeNML(CFG, step, time, v, gamma, psi, alt, lon, lat, cl, cd, mach, temperature, alpha)

########################################
##         Running Propagation        ##
########################################

# Run CFD first to get the coefficients
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
