# Astrophysics: sata pipelines and orbital dynamics

This directory contains Python scripts dedicated to processing astronomical datasets, analyzing stellar spectra and simulating orbital mechanics.

## Key contents

* **Stellar spectrum analysis (`stellar_spectrum_type_*.py`):** computational tools for data cleaning, interpolation and comparison of different stellar spectral types (A, B and K) utilizing specialized scientific libraries.
* **Hertzsprung-Russell diagram analysis (`hr_diagram_analysis_*.py`):** automated pipeline to process stellar datasets, filter experimental noise and structure data to analyze H-R diagrams for Main Sequence (MS) and Supergiant (SG) stars.
* **Data processing pipelines (`astro_data_pipeline_step*.py`):** a multi-step modular pipeline designed to handle, clean and reduce raw experimental data from astronomical observations.
* **Satellite dynamics (`satellite_orbit_dynamics.py`):** numerical simulation modeling the orbital trajectory and stability of a satellite under planetary gravitational constraints.

## Technical stack
* **Python 3**.
* **Astropy:** FITS file handling and astronomical data structuring.
* **SciPy:** numerical interpolation and scientific computing.
* **NumPy & Matplotlib:** matrix manipulation and data visualization.
