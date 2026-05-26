# Electromagnetism: numerical computing and field simulations

This directory contains Python implementations focused on solving complex differential equations, simulating electromagnetic fields and analyzing transient behaviors in dynamic systems.

## Key contents

* **Cyclotron simulation (`cyclotron_magnetic_simulation.py`):** dynamic simulation of charged particles interacting with magnetic fields. Utilizes advanced ODE (Ordinary Differential Equation) solvers to model non-linear trajectories.
* **RLC circuit analysis (`rlc_circuit_transient_analysis_*.py`):** time-series computational model evaluating the transient and steady-state responses of RLC circuits (with and without sources). Highly applicable to signal processing and temporal data modeling.
* **Poisson equation solver (`poisson_equation_solver.py`):** numerical method implementation to solve partial differential equations (PDEs) over spatial grids, calculating and visualizing electrostatic potential distributions.

## Technical stack
* **Python 3**.
* **SciPy:** advanced mathematical integration and ODE solvers (`scipy.integrate`).
* **NumPy:** vectorized operations, matrix manipulations and spatial mesh generation.
* **Matplotlib:** 2D and 3D visualization of physical fields and trajectories.
