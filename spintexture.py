import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import pandas as pd

# reads spin(?) output file
def read_state(filename):
    state = pd.read_csv(filename, sep='\s{2,}', header=None, engine='python')
    k_points = np.array(state.loc[:, '1':'3'])
    energy = state.loc[:, '5':'5']
    energy = np.array(energy).flatten()
    spins = np.array(state.loc[:, '6':'8']).T
    return k_points, energy, spins


def plot_2D_spin_texture(k_points_1_read, energy_1_read, spins_1_read, energy_shift):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams.update({'font.size': 15})
    matplotlib.rc('text', usetex='true')
    fig, ax = plt.subplots(figsize=(8, 5))

    range_set = [i for i in range(441)]

    # normalize the color based on the z-component of the spin
    norm = colors.Normalize(vmin=min(spins_1_read.T[:, 2]), vmax=max(spins_1_read.T[:, 2]))

    for i in range_set:
        color = plt.cm.coolwarm(norm(spins_1_read.T[i][2]))
        q = ax.quiver(k_points_1_read[i][0], k_points_1_read[i][1], spins_1_read.T[i][0], spins_1_read.T[i][1],
                      color=color, zorder=1)

    plt.axvline(x=0, linestyle="--", color="k", linewidth=0.8)
    plt.axhline(y=0, linestyle="--", color="k", linewidth=0.8)

    # create a colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation='vertical', label='$<\sigma_Z>$ component')

    x = k_points_1_read.T[0]
    y = k_points_1_read.T[1]
    z = energy_1_read - energy_shift
    print(min(z))

    # Create a contour plot with specified colormap
    cmap_type = 'gist_earth'  # or 'viridis', 'plasma', 'inferno', 'magma', 'cividis', etc.
    # contour_plot = plt.tricontour(x, y, z, levels=5, vmin=0, vmax=0.1, cmap=cmap_type,
    #                               linestyles='dashed', linewidths=2.0, zorder=2)

    # Add labels to the contour lines
    # plt.clabel(contour_plot, inline=True, fontsize=15)

    # Create a color bar with a Normalize instance
    # cbar = plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=0.1), cmap=cmap_type), ax=plt.gca())
    # cbar.set_label('Contour line - Energy (eV)', fontsize=25)

    # Set the color bar ticks
    # boundaries = np.linspace(0, 0.1, 6)
    # cbar.set_ticks(boundaries)

    plt.xlabel("$k_x$ ($\AA^{-1}$)")
    plt.ylabel("$k_y$ ($\AA^{-1}$)")

    xrange = 0.45
    yrange = 0.23

    plt.xlim(-xrange, xrange)
    plt.ylim(-yrange, yrange)
    plt.title("IJAYUQ VBM (1708)")
    plt.show()


def extra_states(full_file, write_file, state):
    res = []
    with open(full_file, "r") as f:
        for ln in f:
            temp = ln.split()[4]
            if str(state) == temp:
                res.append(ln)
    with open(write_file, "w") as f:
        f.writelines(res)


base = "../../FHI-aims/Double_Perovskites/real-systems/spin_textures/"
# selwoz = base + "selwoz/"
# state = 750
# state = 751
# output = selwoz + "state" + str(state) + ".dat"
# extra_states(selwoz + "spin_texture.dat", output, state)

# selwuf = base + "selwuf/"
# state = 750
# state = 751
# output = selwuf + "state" + str(state) + ".dat"
# extra_states(selwuf + "spin_texture.dat", output, state)

ijayuq = base + "ijayuq/"
state = 1708
# state = 1709
output = ijayuq + "state" + str(state) + ".dat"
extra_states(ijayuq + "spin_texture.dat", output, state)

k_points_1_read, energy_1_read, spins_1_read = read_state(output)
energy_shift = min(energy_1_read)
plot_2D_spin_texture(k_points_1_read, energy_1_read, spins_1_read, energy_shift)
