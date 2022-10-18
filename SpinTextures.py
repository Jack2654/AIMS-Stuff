import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpolate
import math
from os import listdir
from os.path import isfile, join
import math
import seaborn as sns
import pandas as pd
from matplotlib import cm

import matplotlib.ticker as plticker
from mpl_toolkits.mplot3d import Axes3D


def read_state(filename):
    state = pd.read_csv(filename, sep='\s{2,}', header=None, engine='python')
    # print(state)
    k_points = np.array(state.loc[:, '1':'3'])
    energy = state.loc[:, '5':'5']
    energy = np.array(energy).flatten()
    spins = np.array(state.loc[:, '6':'8']).T
    return k_points, energy, spins


k_points_1, energy_1, spins_1 = read_state("/Users/naidelcaturello/Downloads/state_1720_2_spin_texture.dat")
k_points_2, energy_2, spins_2 = read_state("/Users/naidelcaturello/Downloads/state_1721_2_spin_texture.dat")

fig, axs = plt.subplots(3, 3)
fig.set_size_inches(8, 8)
fig.tight_layout(pad=1)

cmpColor = 'bwr'

x_values = [i * 0.3784 / 40 for i in range(48)]  # number of points of k-path
print(len(energy_1), len(k_points_1), len(spins_1))

column = -1
# for segment in [21,42, 63]:
# for segment in [0, 21, 42]:
#    column += 1
for i in range(3):
    pcm = axs[0, i].scatter(x_values, energy_1 - min(energy_1), \
                            c=spins_1[i], cmap=cmpColor, vmin=-0.8, vmax=0.8)

    axs[0, i].scatter(x_values, energy_2 - min(energy_1), \
                      c=spins_2[i], cmap=cmpColor, vmin=-0.8, vmax=0.8)
plt.subplots_adjust(hspace=0.2, wspace=0.25)
# g.set_label(label="Inorganic Contribution", size=12)
plt.show()