# Coupled Trajectory–CFD Controller

**Author:** Kevin Tang  

This repository implements a coupled trajectory–CFD framework for atmospheric entry simulations.  
A Python controller orchestrates trajectory propagation, CFD execution, and state exchange through a shared `config.nml` file, while recording the trajectory history to `history.csv`.

---

## Repository Structure

- `control.py`  
  Python controller coordinating trajectory propagation and CFD execution.

- `Reentry_3DOF_NonPlanar/`  
  Fortran-based 3-DOF non-planar trajectory solver.

- `fun3D_Solver/`  
  FUN3D execution scripts and supporting files.

- `config.nml`  
  Shared namelist file for state exchange between Python and Fortran.

- `history.csv`  
  Time history of trajectory and aerodynamic states.

---

## Instructions

1. Check that all file and directory paths in `control.py` are correct.  
2. Adjust input parameters in `config.nml` as needed.  
3. Ensure `T_STEP` in `config.nml` matches the time stepping used in the Fortran trajectory code.  
4. Enter the ESP shell environment before running the controller.
5. Update the characteristic radius for blunt cone in run_075_70deg.py
6. Activate the Python virtual environment (`.venv`) in the project directory.

---

## GitHub Clone Tip (Library Paths)

If you are using locally built libraries (e.g., NRLMSISE-00), add the following to `~/.bashrc`
or to a project-specific `env.sh` file that can be sourced:

```bash
export CPATH="$HOME/local/include:$CPATH"
export LIBRARY_PATH="$HOME/local/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"
```

Reload the environment after editing:
```bash
source ~/.bashrc
```
