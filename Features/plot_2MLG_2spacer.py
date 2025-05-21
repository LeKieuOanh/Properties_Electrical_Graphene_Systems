import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters
# e22 = [4, 7.5, 22]
e22 = [7.5, 12.53, 22]
d = [3, 6]
w = [1, 2]
# w = [2, 4]

e22 = list(map(str, e22))
d = list(map(str, d))

# Linestyles for r0 values
sb = ['solid', '-.', '--']

# Base path
base_dir = os.getcwd()

# Loop over e22 and d
for j in range(len(e22)):
    for i in range(len(d)):
        # Construct file path
        data_path = os.path.join(
            base_dir,
            'results/MLG_systems_0K/MLG_MLG/2spacer/data',
            f'ni1=0.45e12-d={d[i]}-w={w[i]}-e21=7.5-e22={e22[j]}.txt'
        )

        # Load data
        data = np.loadtxt(data_path, unpack=True)
        x, y1, y2, y3 = data

        # Create figure
        plt.figure(figsize=(6, 5))
        plt.xlabel('n1', fontsize=10)
        plt.ylabel('$\sigma(e^2/h)$', fontsize=15)
        plt.title(f'2MLG - ni1=0.45e12 - d = {d[i]}nm - w = {w[i]}nm - e21=7.5 - e22 = {e22[j]}')

        # Plot with distinct styles
        plt.plot(x, y1, ls=sb[0], label='r0 = 0', color='blue')
        plt.plot(x, y2, ls=sb[1], label='r0 = 5a0', color='green')
        plt.plot(x, y3, ls=sb[2], label='r0 = 10a0', color='red')

        plt.legend(fontsize=8)
        plt.tight_layout()

        # Save figure
        save_path = os.path.join(
            base_dir,
            'results/MLG_systems_0K/MLG_MLG/2spacer/image',
            f'd={d[i]}-w={w[i]}-ni1=0.45e12-e21=7.5-e22={e22[j]}.pdf'
        )
        plt.savefig(save_path, dpi=1000)
        plt.close()  # Avoid memory overload when looping
