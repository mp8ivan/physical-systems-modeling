# Quantum mechanics: linear algebra and sparse matrices

This directory highlights computational approaches to solving complex quantum mechanical problems, emphasizing the use of advanced linear algebra, matrix operations and numerical comparisons.

## Key contents

* **Quantum state evolution (`quantum_state_evolution_1d_*.py`):** algorithmic simulation modeling the time evolution of different initial wave functions with a harmonic potential and with an infinite potential well. Implements sparse matrix representations to optimize memory and processing speed for large-scale calculations.
* **Analytical vs numerical models (`numerical_vs_analytical_quantum_*.py`):** a robust comparison script validating numerical approximation methods against exact analytical solutions, demonstrating rigorous mathematical grounding and error analysis.

## Technical stack
* **Python 3**.
* **SciPy (`scipy.sparse`):** efficient implementation and manipulation of large sparse matrices, crucial for optimizing computational resources in high-dimensional state spaces.
* **NumPy:** core linear algebra operations and eigenvalue/eigenvector computations.
* **Matplotlib:** visualization of wave functions and probability density distributions.
