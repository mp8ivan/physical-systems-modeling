# Electromagnetism: numerical computing and field simulations

This directory contains Python implementations focused on solving complex differential equations, simulating electromagnetic fields and analyzing transient behaviors in dynamic systems.

## Key contents

* **Cyclotron simulation (`cyclotron_magnetic_simulation.py`):** dynamic simulation of charged particles interacting with magnetic fields. Utilizes advanced ODE (Ordinary Differential Equation) solvers to model non-linear trajectories.
* **RLC circuit analysis (`rlc_circuit_transient_analysis_*.py`):** time-series computational model evaluating the transient and steady-state responses of RLC circuits (with and without sources). Highly applicable to signal processing and temporal data modeling.
* **Poisson equation solver (`poisson_equation_solver.py`):** numerical method implementation to solve partial differential equations (PDEs) over spatial grids, calculating and visualizing electrostatic potential distributions.
* **Square charge poisson solver (`square_charge_poisson_solver.py`):** a specialized electrostatic simulation that evaluates potential fields under a discrete square charge distribution or boundary condition. Demonstrates proficiency in applying relaxation methods and matrix algebra to tailored geometric constraints.

## Technical stack
* **Python 3**.
* **SciPy:** advanced mathematical integration and ODE solvers (`scipy.integrate`).
* **NumPy:** vectorized operations, matrix manipulations and spatial mesh generation.
* **Matplotlib:** 2D and 3D visualization of physical fields and trajectories.
