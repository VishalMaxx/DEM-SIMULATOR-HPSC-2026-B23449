import matplotlib.pyplot as plt
import numpy as np

# Load the data from your simulator output
data = np.loadtxt('results.txt', skiprows=1)
time = data[:, 0]
num_z = data[:, 1]
analytic_z = data[:, 2]
ke = data[:, 3]

# --- Plot 1: Trajectory Comparison (Verification) ---
plt.figure(figsize=(8, 5))
plt.plot(time, num_z, 'b-', label='Numerical (DEM)', linewidth=2)
plt.plot(time, analytic_z, 'r--', label='Analytical', linewidth=2)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Height Z (m)', fontsize=12)
plt.title('Test 1: Free Fall Verification', fontsize=14)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('trajectory_plot.png', dpi=300)
print("Saved trajectory_plot.png")

# --- Plot 2: Kinetic Energy (Damping Check) ---
plt.figure(figsize=(8, 5))
plt.plot(time, ke, 'g-', linewidth=2)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Kinetic Energy (J)', fontsize=12)
plt.title('Kinetic Energy vs Time (Bounce Test)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('energy_plot.png', dpi=300)
print("Saved energy_plot.png")
