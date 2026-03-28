# Parallel DEM Simulator (HPSC Assignment 1)
**Author:** Vishal (IIT Mandi)

## Features
- **Physics:** 3D Discrete Element Method with Spring-Dashpot contact laws.
- **Parallelization:** OpenMP multi-threading (Scales up to 12+ threads).
- **Optimization:** $O(N)$ Neighbor Search using Cell-Linked Lists.
- **Verification:** Analytical comparison for Free Fall and Kinetic Energy decay.


## Simulation Results

### 1. Free Fall Verification
Comparison between numerical DEM results and the analytical solution.
![Trajectory](results/plots/trajectory_plot.png)

### 2. Energy Dissipation
Demonstration of the Spring-Dashpot damping effect.
![Energy](results/plots/energy_plot.png)

### 3. Scientific Studies
Timestep sensitivity (Convergence) and Damping effect analysis.
![Sensitivity](results/plots/sensitivity_plot.png)
![Damping](results/plots/damping_plot.png)

### 4. Strong Scaling Analysis
Performance scaling from 1 to 12 threads using OpenMP, demonstrating the impact of Amdahl's Law on the $O(N)$ optimized engine.
![Scaling](results/plots/scaling_plot.png)


## Repository Branch Structure
To clearly separate the baseline requirements from the bonus optimization tasks, this repository is organized into three branches:

* **`main`**: Contains the final, fully integrated codebase (combining all optimizations and studies).
* **`neighbor-search`**: Isolates the commits where the $O(N)$ Cell-Linked List spatial partitioning was implemented (Section 18).
* **`scientific-study`**: Isolates the commits containing the timestep sensitivity and damping coefficient scripts (Section 19).
