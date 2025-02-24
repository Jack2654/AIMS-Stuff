import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import pandas as pd
import BasicBandOut as bbo


# reads stateXXXX.out output file
def read_state(filename):
    state = pd.read_csv(filename, sep='\s{2,}', header=None, engine='python')
    k_points = np.array(state.loc[:, '1':'3'])
    energy = state.loc[:, '5':'5']
    energy = np.array(energy).flatten()
    spins = np.array(state.loc[:, '6':'8']).T
    return k_points, energy, spins


def plot_2D_spin_texture(k_points_1_read, energy_1_read, spins_1_read, energy_shift, filename, title="", save=False):
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
    fig.colorbar(sm, ax=ax, orientation='vertical', label=r'$\sigma_{\mbox{z}}$ component')

    # x = k_points_1_read.T[0]
    # y = k_points_1_read.T[1]
    # z = energy_1_read - energy_shift
    # Create a contour plot with specified colormap
    # cmap_type = 'gist_earth'  # or 'viridis', 'plasma', 'inferno', 'magma', 'cividis', etc.
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

    plt.xlabel(r'$k_x$ (\AA$^{-1}$)')
    plt.ylabel(r'$k_y$ (\AA$^{-1}$)')

    xrange = 0.45
    yrange = 0.45

    plt.xlim(-xrange, xrange)
    plt.ylim(-yrange, yrange)
    plt.title(title)
    if save:
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()
    plt.clf()


def extra_states(full_file, write_file, state):
    res = []
    with open(full_file, "r") as f:
        for ln in f:
            temp = ln.split()[4]
            if str(state) == temp:
                res.append(ln)
    with open(write_file, "w") as f:
        f.writelines(res)


base = "../../FHI-aims/Double_Perovskites/Real_Structures/real-systems/spin_textures/"
structure = base + "selwoz/"
states = [749, 750, 751, 752]  # VBM lower and upper, CBM lower and upper

# structure = base + "selwuf/"
# states = [749, 750, 751, 752]  # VBM lower and upper, CBM lower and upper

# structure = base + "ijayuq/" # limits 0.40, 0.225
# states = [1707, 1708, 1709, 1710]  # VBM lower and upper, CBM lower and upper

for state in states:
    # continue
    cur_state = structure + "state" + str(state)
    output = cur_state + ".dat"
    extra_states(structure + "spin_texture.dat", output, state)

    k_points_1_read, energy_1_read, spins_1_read = read_state(output)
    energy_shift = min(energy_1_read)
    title = "SELWOZ "
    if state % 2 == 0:
        title += "Upper "
    else:
        title += "Lower "
    if state < sum(states)/4:
        title += "VBM "
    else:
        title += "CBM "
    title += f'(%s)' % state

    plot_2D_spin_texture(k_points_1_read, energy_1_read, spins_1_read, energy_shift, cur_state + ".png", title=title,
                         save=False)

base = "../../FHI-aims/Double_Perovskites/Real_Structures/real-systems/"
options = ['experimental-bs-selwoz/', 'experimental-bs-selwuf/', 'ijayuq/']
for opt in options:
    print(opt)
    for i in range(1, 2):
        print(bbo.band_info(base + opt, f'band100{i}.out', band_gap=True, spin_splitting=True, verbose=False))
