# Classical mechanics: chaotic systems and non-linear dynamics

This directory focuses on the computational modeling of classical physical systems, highlighting the numerical resolution of non-linear dynamics, rigid body rotations and non-inertial reference frames.

## Key contents

* **Double pendulum chaos (`double_pendulum_chaos.py`):** an advanced numerical simulation of a double pendulum system solving coupled non-linear ordinary differential equations (ODEs) to model highly sensitive, chaotic motion.
* **Damped double pendulum (`damped_double_pendulum.py`):** an extension of the chaotic pendulum model that incorporates physical damping and friction factors, demonstrating energy dissipation in dynamic systems.
* **Tennis racket theorem (`tennis_racket_theorem.py`):** computational modeling of the Euler equations for rigid body dynamics. It simulates the intermediate axis theorem (Dzhanibekov effect), calculating moments of inertia and precession frequencies in 3D space.
* **Coriolis projectile motion (`coriolis_projectile_motion.py`):**a 3D kinematics simulation of a projectile that accounts for aerodynamic drag and the Coriolis effect by factoring in the Earth's rotation, latitude and firing azimuth.

## Technical stack
* **Python 3**.
* **SciPy:** advanced ODE integration (`scipy.integrate.odeint`) for solving coupled differential equations and rigid body rotations.
* **NumPy:** state-vector manipulation, array broadcasting and trigonometric transformations.
* **Matplotlib:** kinematic animation, phase-space trajectory visualization and 3D plotting (`mpl_toolkits.mplot3d`).
