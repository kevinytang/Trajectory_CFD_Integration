# Import pyCAPS module
import pyCAPS
import os
import numpy as np
import f90nml
import shutil


def fbool(x: bool) -> str:
    return ".true." if bool(x) else ".false."


#######################################
##        Load input parameters      ##
#######################################
params = f90nml.read("../config.nml")
cfdvar = f90nml.read("cfd_control.nml")

# Add access to each sections
traj = params["trajectory_input"]
cfd  = params["cfd_input"]
var  = cfdvar["cfd_var"]

#######################################
##        Build Geometry             ##
#######################################
filename = "075_70deg.csm"
print(f'\n==> Loading geometry from file "{filename}"...')
capsProblem = pyCAPS.Problem(problemName="075_70deg",
                             capsFile=filename,
                             outLevel=1)

#######################################
##        Build surface mesh         ##
#######################################
print('\n==> Creating AFLR4 AIM')
aflr4 = capsProblem.analysis.create(aim="aflr4AIM", name="aflr4")

# Aflr4 inputs
aflr4.input.Mesh_Length_Factor = var.get("Mesh_Length_Factor")
aflr4.input.max_scale = var.get("max_scale")
aflr4.input.ideal_min_scale = var.get("min_scale")
aflr4.input.ff_cdfr = var.get("ff_cdfr")

# Assign geometric groups - for inviscid, all surfaces are inviscid
aflr4.input.Mesh_Sizing = {
    "blunt": {"edgeWeight": var.get("edgeWeight"),
              "scaleFactor": var.get("blunt_scaleFactor")},
    "Farfield": {
        "bcType": "Farfield",
        "scaleFactor": var.get("farfield_scaleFactor")
    }
}

print('\n==> Running AFLR4 (pre/post-analysis)')
aflr4.runAnalysis()

#######################################
##        Build volume mesh          ##
#######################################
print('\n==> Creating AFLR3 AIM')
aflr3 = capsProblem.analysis.create(aim="aflr3AIM", name="aflr3")

# Link AFLR4 surface mesh to AFLR3 (parent/child)
aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

# Define groups: mark walls as Inviscid (no BL layers)
aflr3.input.Mesh_Sizing = {
    "blunt": {"bcType": "Inviscid"},
    "Farfield": {"bcType": "Farfield"}
}

print("==> Running AFLR3 (pre/post-analysis)")
aflr3.runAnalysis()

#######################################
##         Using FUN3D AIM          ##
#######################################
print('\n==> Creating FUN3D AIM')
fun3d = capsProblem.analysis.create(aim="fun3dAIM", name="fun3d")

# Project name & mesh link
fun3d.input.Proj_Name = "075_70deg"
fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])

# Namelist generation from Python
fun3d.input.Use_Python_NML = var.get("Use_Python_NML", True)
fun3d.input.Overwrite_NML = var.get("Overwrite_NML", True)

# MPI procs
nproc = int(os.environ.get("FUN3D_MPI_PROCS", str(var.get("np"))))

# Iterations and CFL schedule
fun3d.input.Num_Iter = int(var.get("Num_Iter"))
fun3d.input.CFL_Schedule = var.get("CFL_Schedule")
fun3d.input.CFL_Schedule_Iter = var.get("CFL_Schedule_Iter")
fun3d.input.Restart_Read = var.get("Restart_Read")

# Boundary conditions for inviscid hypersonic flow
fun3d.input.Boundary_Condition = {
    "blunt": {
        "bcType": "Inviscid"
    },
    "Farfield": {
        "bcType": "Freestream",
        "machNumber": cfd.get("mach")
    }
}

########## Run FUN3D ##########
print("\n\nRunning FUN3D......")

# Write AIM-generated namelist first
fun3d.preAnalysis()

# Compute characteristic surface area for spherically blunted cone
r = 2  # m
A_total = np.pi * r ** 2

# Store the surface length
var["Characteristic_Surface"] = float(A_total)

# ---- Patch extra blocks into fun3d.nml ----
nml_path = os.path.join(fun3d.analysisDir, "fun3d.nml")
with open(nml_path, "a") as f:
    # ---- Reference_physical_properties ----
    f.write("&reference_physical_properties\n")
    f.write(f"  dim_input_type = '{var.get('dim_input_type')}'\n")
    f.write(f"  temperature_units = '{var.get('temperature_units')}'\n")
    f.write(f"  gridlength_conversion = {float(var.get('gridlength_conversion'))}\n")
    f.write(f"  mach_number = {float(cfd.get('mach'))}\n")
    f.write(f"  temperature = {float(cfd.get('temperature'))}\n")
    f.write(f"  angle_of_attack = {float(cfd.get('alpha'))}\n")
    f.write(f"  angle_of_yaw = {float(var.get('Beta'))}\n")
    f.write("/\n\n")

    # ---- Governing equations ----
    f.write("&governing_equations\n")
    f.write(f"  eqn_type = '{var.get('Equation_Type')}'\n")
    f.write(f"  viscous_terms = '{var.get('Viscous')}'\n")
    f.write("/\n\n")

    # ---- Inviscid flux method (CRITICAL for hypersonic) ----
    f.write("&inviscid_flux_method\n")
    f.write(f"  flux_construction = '{var['Flux_Construction']}'\n")
    f.write(f"  flux_limiter = '{var['Flux_Limiter']}'\n")
    f.write(f"  first_order_iterations = {int(var['First_Order_Iterations'])}\n")
    f.write(f"  freeze_limiter_iteration = {int(var['Freeze_Limiter_Iteration'])}\n")
    f.write("/\n\n")

    # ---- Nonlinear solver parameters ----
    f.write("&nonlinear_solver_parameters\n")
    f.write(f"  schedule_iteration = {var['CFL_Schedule_Iter'][0]} {var['CFL_Schedule_Iter'][1]}\n")
    f.write(f"  schedule_cfl = {var['CFL_Schedule'][0]} {var['CFL_Schedule'][1]}\n")
    f.write("/\n\n")

    # ---- Reference Area and Lengths ----
    f.write("&force_moment_integ_properties\n")
    f.write(f"  area_reference = {var['Characteristic_Surface']}\n")
    f.write("/\n\n")

    # Linear solver control block
    f.write("&linear_solver_control\n")
    f.write("  linear_projection = .true.\n")
    f.write("/\n\n")

# Path to nodet_mpi
nodet = var.get('nodet_Dir')
cmd = f"mpirun -np {nproc} {nodet} --animation_freq -1 --volume_animation_freq -1"

# Freeze limiter (optional but recommended for hypersonic)
freeze_limiter = var.get("Freeze_Limiter")
if freeze_limiter:
    cmd += f" --freeze_limiter {freeze_limiter}"

print("Command:", cmd)
fun3d.system(cmd)

# Post-analysis
fun3d.postAnalysis()
print("\nDone. Check Info.out, fun3d.nml, and mapbc.dat.")

#######################################
##         Extract CL and CD         ##
#######################################
force_file = os.path.join(fun3d.analysisDir, f"{fun3d.input.Proj_Name}.forces")

# Pre-allocate
CL, CD = None, None

# Extract CL and CD
if os.path.exists(force_file):
    with open(force_file) as f:
        lines = f.readlines() # read all lines into a list

        # Find last Cl and Cd appearance
        for line in reversed(lines):
            if line.strip().startswith("Cl") and "Cd" in line:

                # Clean up the line
                line = line.strip()

                # Replace = with space
                parts = line.replace("=", " ").split()

                # Extract numerical values
                CL = float(parts[1])
                CD = float(parts[3])
                break

# Update the namelist
if CL is not None and CD is not None:
    traj["cl"] = CL
    traj["cd"] = CD
    with open("../config.nml", "w") as f:
        f90nml.write(params, f)
    print("Updated ../config.nml with CL/CD.")
else:
    print("WARNING: Could not extract CL/CD from .forces file.")