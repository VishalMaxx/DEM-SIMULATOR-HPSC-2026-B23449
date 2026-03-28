import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt('damping_results.txt')
gammas = np.unique(data[:, 0])

plt.figure(figsize=(10, 6))

for g in gammas:
    subset = data[data[:, 0] == g]
    plt.plot(subset[:, 1], subset[:, 2], label=f'Gamma = {g}')

plt.axhline(y=0.5, color='k', linestyle='--', alpha=0.3, label='Floor (Radius)')
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Height Z (m)', fontsize=12)
plt.title('Section 19.1: Effect of Damping on Rebound Height', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('damping_plot.png', dpi=300)
print("Saved damping_plot.png")
