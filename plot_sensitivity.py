import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt('sensitivity_results.txt')
dts = data[:, 0]
errors = data[:, 1]

plt.figure(figsize=(8, 6))
# Log-Log plot is best for convergence studies
plt.loglog(dts, errors, 'bo-', label='Measured Error', markersize=8)

# Add a reference line for O(dt) convergence
plt.loglog(dts, dts * (errors[0]/dts[0]), 'k--', label='Order 1 Reference')

plt.xlabel('Timestep size (dt)', fontsize=12)
plt.ylabel('Absolute Error at t=1.0s (m)', fontsize=12)
plt.title('Timestep Sensitivity Study (Convergence)', fontsize=14)
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()
plt.savefig('sensitivity_plot.png', dpi=300)
print("Saved sensitivity_plot.png")
