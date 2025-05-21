import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters
e22 = [12.53, 22]
d = [3, 6]
e22 = list(map(str, e22))
d = list(map(str, d))

linestyles = ['solid', '-.']  # Different for e22[0] and e22[1]
colors = ['blue', 'green', 'red']  # For y1, y2, y3

# Base path
base_dir = os.getcwd()

# Create figure with 2 rows (for d) and 2 columns (for e22)
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharex=True, sharey=True)

# Loop and plot each file in a subplot
for j in range(len(e22)):  # e22 → columns
    for i in range(len(d)):  # d → rows
        # Load data
        data_path = os.path.join(
            base_dir,
            'results/MLG_systems_0K/MLG_MLG/2spacer/data',
            f'd={d[i]}-ni1=0.45e12-e21=7.5-e22={e22[j]}.txt'
        )
        data = np.loadtxt(data_path, unpack=True)
        x, y1, y2, y3 = data

        # Get subplot axes
        ax_ij = ax[i, j]
        ax_ij.plot(x, y1, color=colors[0], label='r0 = 0')
        ax_ij.plot(x, y2, color=colors[1], label='r0 = 5a0')
        ax_ij.plot(x, y3, color=colors[2], label='r0 = 10a0')

        # Title and labels
        ax_ij.set_title(f'd={d[i]}nm, e22={e22[j]}', fontsize=10)
        ax_ij.set_xlabel('n1')
        ax_ij.set_ylabel('$\sigma(e^2/h)$')
        ax_ij.legend(fontsize=7)

# Tight layout
plt.tight_layout()

# Save combined figure
save_path = os.path.join(
    base_dir,
    'results/MLG_systems_0K/MLG_MLG/2spacer/image',
    'combined_e21=7.5_2MLG_subplots.pdf'
)
plt.savefig(save_path, dpi=600)
plt.show()
