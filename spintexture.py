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
    plt.tight_layout()
    if save:
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()
    plt.clf()


# plot_2D_spin_texture(k_points_1_read, energy_1_read, spins_1_read, energy_shift, cur_state + ".png", title=title, save=False)


def extra_states(full_file, write_file, state):
    res = []
    with open(full_file, "r") as f:
        for ln in f:
            temp = ln.split()[4]
            if str(state) == temp:
                res.append(ln)
    with open(write_file, "w") as f:
        f.writelines(res)


def plot_all(x_loc, y_loc, x_spin, y_spin, z_spin):
    nrows, ncols = 2, 4
    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(16, 8))
    axes = axes.flatten()
    norm = colors.Normalize(vmin=min([min(sv) for sv in z_spin]), vmax=max([max(sv) for sv in z_spin]))

    for ax, x, y, a, b, c in zip(axes, x_loc, y_loc, x_spin, y_spin, z_spin):
        for i in range(len(x)):
            color = plt.cm.coolwarm(norm(c[i]))
            ax.quiver(x[i], y[i], a[i], b[i], color=color, zorder=1)
        ax.set_xlabel(r"$k_x$")
        ax.axvline(x=0, linestyle="--", color="k", linewidth=0.8)
        ax.axhline(y=0, linestyle="--", color="k", linewidth=0.8)
    axes[0].set_ylabel(r"$k_y$")
    axes[4].set_ylabel(r"$k_y$")
    titles = [r"(OCA)$_4$AgBiBr$_8$", r"(R-3AP)$_4$AgBiBr$_8$", r"(S-3AP)$_4$AgBiBr$_8$",
              r"(R-$\beta$-MPA)$_4$AgBiI$_8$"]

    for col, title in enumerate(titles):
        axes[col].set_title(title, pad=10)

    for title, y in zip(["Lower CBM", "Upper CBM"], [0.28, 0.68]):
        fig.text(0.02, y, title, va="center", ha="center", fontsize=12, rotation="vertical")

    for i, ax in enumerate(axes):
        col = i % ncols
        if col != 0:
            ax.tick_params(left=False, labelleft=False)
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.08, top=0.88, wspace=0.0, hspace=0.0)
    cax = fig.add_axes([0.90, 0.08, 0.02, 0.80])
    sc = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=norm)
    sc.set_array([])
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label(r"$\sigma_z$")
    plt.show()


base = "../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/spin_textures/"
structures = ["voxkif", "selwoz", "selwuf", "ijayuq"]
vbms = [[709, 710], [749, 750], [749, 750], [1707, 1708]]
cbms = [[712, 711], [752, 751], [752, 751], [1710, 1709]]
x_loc, y_loc, x_spin, y_spin, z_spin = [], [], [], [], []
for i in range(2):
    for j in range(4):
        cur_base = base + structures[j] + "/"
        state = cbms[j][i]
        cur_state = cur_base + "state" + str(state)
        output = cur_state + ".dat"
        extra_states(cur_base + "spin_texture.dat", output, state)

        k_points_1_read, energy_1_read, spins_1_read = read_state(output)
        x_loc.append(k_points_1_read[:, 0])
        y_loc.append(k_points_1_read[:, 1])
        x_spin.append(spins_1_read.T[:, 0])
        y_spin.append(spins_1_read.T[:, 1])
        z_spin.append(spins_1_read.T[:, 2])

plot_all(x_loc, y_loc, x_spin, y_spin, z_spin)
