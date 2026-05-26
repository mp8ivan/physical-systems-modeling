# Quantum mechanics and solid state: linear algebra and lattice structures

This directory highlights computational approaches to solving complex quantum mechanical and structural solid-state problems, emphasizing the use of advanced linear algebra, matrix operations and spatial mesh generation.

## Key contents

* **Quantum state evolution (`quantum_state_evolution_1d_*.py`):** algorithmic simulation modeling the time evolution of different initial wave functions. Implements sparse matrix representations to optimize memory and processing speed for large-scale calculations.
* **Analytical vs numerical models (`numerical_vs_analytical_quantum_*.py`):** a robust comparison script validating numerical approximation methods against exact analytical solutions, demonstrating rigorous mathematical grounding and error analysis.
* **Diamond cubic lattice with basis (`diamond_cubic_lattice_basis.py`):** a crystallography simulation that models and visualizes a diamond cubic crystal structure using conventional lattice vectors and an atomic basis. It bridges quantum mechanics and solid-state physics by demonstrating spatial coordinate transformations and periodic boundary representations in 3D space.

## Technical stack
* **Python 3**.
* **SciPy (`scipy.sparse`):** efficient implementation and manipulation of large sparse matrices, crucial for optimizing computational resources in high-dimensional state spaces.
* **NumPy:** core linear algebra operations, matrix manipulations and 3D coordinate transformations.
* **Matplotlib:** visualization of wave functions, probability densities and 3D crystal lattice projections (`mpl_toolkits.mplot3d`).
