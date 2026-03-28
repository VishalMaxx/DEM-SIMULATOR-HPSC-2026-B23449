import matplotlib.pyplot as plt

# Data from your terminal runs (N=5000)
threads = [1, 2, 4, 8, 12]
times = [6.1476, 4.6820, 2.9974, 2.3196, 1.8422]

# Calculate Speedup (T1 / Tp) and Efficiency (Speedup / p)
speedup = [times[0] / t for t in times]
efficiency = [(s / p) * 100 for s, p in zip(speedup, threads)]

fig, ax1 = plt.subplots(figsize=(8, 5))

# Plot 1: Speedup (Left Axis)
color = 'tab:blue'
ax1.set_xlabel('Number of Threads', fontsize=12)
ax1.set_ylabel('Speedup (S)', color=color, fontsize=12)
ax1.plot(threads, speedup, 'bo-', linewidth=2, label='Actual Speedup')
ax1.plot(threads, threads, 'k--', alpha=0.5, label='Ideal Speedup') # Ideal line
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xticks(threads)
ax1.legend(loc='upper left')

# Plot 2: Efficiency (Right Axis)
ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Efficiency (%)', color=color, fontsize=12)  
ax2.plot(threads, efficiency, 'rs-', linewidth=2, label='Efficiency')
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0, 110)
ax2.legend(loc='upper right')

plt.title('Strong Scaling Analysis (N=5000)', fontsize=14)
plt.grid(True, alpha=0.3)
plt.savefig('scaling_plot.png', dpi=300, bbox_inches='tight')
print("Saved scaling_plot.png")
