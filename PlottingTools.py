import math

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import BasicGeo as bg
import BasicControl as bc
import BasicAimsOut as bao
import BasicBandOut as bbo
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
import matplotlib.patheffects as pe
import scipy
import Geometry


# reimagined arguments: python3 *.py a1 a2 a3 a4 a5 a6 a7 a8 a9
# a1 = file path to folder containing band files
# a2 = name of figure wanted to be saved
# a3 = energyshift
# a4 = ymin
# a5 = ymax
# a6 = substate
# a7 = color dictionary
# a8 = equal sized plotted band segments
# a9 = debug

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, 'arrow3D', _arrow3D)


def dos_plot(folder, shift=0, limits=None, combine=True, save=False, title="Species Projected DOS", filename=None):
    if limits is None:
        limits = [-10, 10]
    all_files = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    files = []
    for file in all_files:
        if ("_l_proj_dos.dat" in file) and ("no_soc" not in file):
            files.append(folder + file)
    fig = plt.gcf()
    fig.set_size_inches(4, 6)
    matplotlib.rc('text', usetex='true')
    matplotlib.rcParams['axes.linewidth'] = 2
    matplotlib.rcParams.update({'font.size': 20})
    colorkey = {"Pb": "m", "Bi": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
    ax = plt.axes()
    doses = []
    elements = []
    for file in files:
        lines = []
        with open(file, "r") as f:
            for ln in f:
                if "#" not in ln:
                    lines.append([float(x) for x in ln.split()[:2]])
        energy = []
        dos = []
        for line in lines:
            if limits[0] < line[0] + shift < limits[1]:
                energy.append(line[0] + shift)
                dos.append(line[1])
        element = file.split("/")[-1].split("_")[0]
        doses.append(dos)
        elements.append(element)
    ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    plt.ylabel("Energy (eV)")
    plt.title(title)
    if combine:
        total_dos = [0 for x in doses[0]]
        organic_dos = [0 for x in doses[0]]
        for index, element in enumerate(elements):
            if element == "Pb" or element == "I" or element == "Bi":
                ax.plot(doses[index], energy, color=colorkey[element], label=element)
            else:
                organic_dos = [organic_dos[i] + doses[index][i] for i in range(len(organic_dos))]
            total_dos = [total_dos[i] + doses[index][i] for i in range(len(total_dos))]
        ax.plot(organic_dos, energy, color='b', label="Organic", alpha=0.5)
        ax.fill_betweenx(energy, 0, total_dos, color='k', alpha=0.1, label="Total DOS")
    else:
        for index, element in enumerate(elements):
            ax.plot(doses[index], energy, color=colorkey[element], label=element)
    plt.legend(fontsize="15")
    plt.tight_layout()
    if save:
        if filename is None:
            filename = "/".join(files[0].split("/")[:-1]) + "/dos.png"
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
        plt.clf()
    else:
        plt.show()


def MD_plots(file):
    lines = []
    with open(file, "r") as f:
        for count, ln in enumerate(f):
            if "#" not in ln:
                lines.append(ln)
    times = [float(line.split()[0]) for line in lines]
    temps = [float(line.split()[1]) for line in lines]
    energies = [float(line.split()[2]) for line in lines]
    energies = [x - energies[0] for x in energies]

    # temperature over time
    plt.plot(times, temps)
    plt.xlabel("Time (ps)")
    plt.ylabel("Temperature (K)")
    plt.ylim([0, 620])
    plt.show()
    plt.clf()

    # energy over time
    plt.plot(times, energies)
    plt.xlabel("Time (ps)")
    plt.ylabel("Change in Energy (eV)")
    plt.show()


def many_MD(folder):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    legend = []
    for direct in all_dir:
        if os.path.isdir(folder + direct):
            data = folder + direct + "/MD.dat"
            lines = []
            legend.append(direct)
            with open(data, "r") as f:
                for ln in f:
                    if "#" not in ln:
                        lines.append(ln)

            times = [float(line.split()[0]) for line in lines]
            temps = [float(line.split()[1]) for line in lines]
            energies = [float(line.split()[2]) for line in lines]

            plt.plot(times, temps)
            plt.xlabel("Time (ps)")
            plt.ylabel("Temperature (K)")
            # plt.ylabel("Energy (eV)")
            plt.ylim([0, 620])
            plt.legend(legend)
            plt.title("Temperature over time for H5O2 MD")
    plt.show()


# should add check that number of labels matches number of bands upfront
def mulliken_plot_old(settings_file, alt_geo=False, debug=False, quiet=False, save=False):
    ############################################
    # Setup                                    #
    ############################################
    input_settings = read_band_settings(settings_file, debug)
    filepath = input_settings[0]
    filename = input_settings[1]
    energyshift = input_settings[2]
    ymin = input_settings[3]
    ymax = input_settings[4]
    substate = input_settings[5]
    color_dict = input_settings[6]
    labels = input_settings[7]
    title = input_settings[8]
    eq = input_settings[9]
    bands = input_settings[10]
    mod = int(input_settings[11])

    flags = 1  # normal run state, print out timing
    if debug:
        flags = 0  # full debug state, print all info
    elif quiet:
        flags = 2  # quiet run state, print out nothing

    ############################################
    # geometry.in                              #
    ############################################
    if not alt_geo:
        rlatvec, atoms, species, species_id = bg.read_geo_for_bands(filepath, color_dict, flags)
    else:
        rlatvec, atoms, species, species_id = bg.read_geo_for_bands(alt_geo, color_dict, flags)
    ############################################
    # control.in                               #
    ############################################
    xvals, band_len_tot, tot_len = bc.read_control_for_bands(filepath, rlatvec, bands, eq, debug)

    ############################################
    # pre-processing                           #
    ############################################
    species_color = color_dict
    if not species_color:
        python_colors = ['r', 'y', 'b', 'c', 'm', 'g', 'pink', 'purple', 'black']
        for s in range(len(species)):
            species_color[species[s]] = python_colors[s]

    sub_st = False
    if not substate == 0:
        sub_st = True  # 0, s; 1, p; 2, d; ...

    linewidth = 1
    black_bands = 1
    markersizeunit = 2
    line_mlk_start_id = 5  # column of contributions in mulliken output
    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    matplotlib.rc('text', usetex='true')
    matplotlib.rcParams['axes.linewidth'] = 2
    matplotlib.rcParams.update({'font.size': 20})

    all_files = [join(filepath, f) for f in listdir(filepath) if isfile(join(filepath, f))]
    band_files = [f for f in all_files if 'band' in f.split("/")[-1] and f.endswith('.out') and 'bandmlk' not in f]
    band_files.sort()
    band_files = [band_files[x] for x in bands]
    band_mlkfiles = [f for f in all_files if 'bandmlk' in f and f.endswith('.out')]
    band_mlkfiles.sort()
    mlk = False
    if len(band_mlkfiles):
        mlk = True
        band_mlkfiles = [band_mlkfiles[x] for x in bands]

    if not flags == 2:
        print("color key: ", species_color)
        if flags == 0:
            print(band_files)
            if mlk:
                print(band_mlkfiles)

    ############################################
    # band10**.out							   #
    ############################################
    if not flags == 2:
        print("Analyzing band10**.out files...")
    energys = []
    for file_id in range(len(band_files)):
        temp_en = []
        file = band_files[file_id]
        curTime = time.time()
        if not flags == 2:
            print("Processing " + file, end='\t\t')
        f = open(file)
        for line in f:
            words = line.strip().split()
            energy = []
            occ_ene = words[4:]
            for i in range(0, len(occ_ene), 2):
                energy.append(float(occ_ene[i + 1]) - energyshift[file_id])
            temp_en.append(energy)
        # len(energy) = number of states; len(temp_en) = number of k points
        # band_points.append(len(k_index))
        for i in range(len(energy)):
            band = []
            for j in range(len(temp_en)):
                band.append(temp_en[j][i])
            plt.plot(xvals[file_id], band, color='k', lw=black_bands)
        energys.append(temp_en)
        if not flags == 2:
            print(str(time.time() - curTime) + " seconds")

    ############################################
    # output non-mlk						   #
    ############################################
    if not mlk:
        x_pts = []
        x_pts.append(0)
        for kpoint_x in band_len_tot[1:]:
            x_pts.append(kpoint_x)
            plt.axvline(kpoint_x, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
            plt.axhline(0, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])

        x_pts.append(tot_len)
        if labels == 0:
            print("No labels found")
        else:
            plt.xticks(ticks=x_pts, labels=labels)
        plt.ylabel('Energy (eV)')
        plt.xlabel("Wave vector, $k$")
        plt.title(title, weight='bold')
        plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
        plt.tight_layout()
        if save and not filename == 0:
            plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
        else:
            plt.show()
        plt.clf()
        return 0

    ############################################
    # mlk									   #
    ############################################
    if not flags == 2:
        print("Analyzing bandmlk100*.out files...")

    # for file_id in range(len(band_mlkfiles)):
    for file_id, mlkfile in enumerate(band_mlkfiles):
        mlkfile = band_mlkfiles[file_id]
        curTime = time.time()
        print("Processing " + mlkfile, end='\t\t')
        f = open(mlkfile)
        currentMlk = []
        i = 0
        for line in f:
            if line.strip() == '':
                continue
            elif line.strip().startswith("k point number:"):
                currentK = int(line.strip().split()[3][:-1])
                kpt = line.strip()
                currentMlk = []
            elif line.strip().startswith("State") and len(line.strip().split()) == 2:
                words = line.strip().split()
                currentState = int(words[1])
                currentMlk = []
            elif line.strip()[0].isdigit():
                words = line.strip().split()
                if sub_st:
                    currentMlk.append(float(words[line_mlk_start_id + 1 + substate]))
                else:
                    only_pos = False
                    if only_pos:
                        temp = 0
                        for x in words[line_mlk_start_id + 1:]:
                            if float(x) > 0:
                                temp += float(x)
                        currentMlk.append(temp)
                    else:
                        currentMlk.append(float(words[line_mlk_start_id]))
            current_energy = float(words[1])
            if len(currentMlk) == len(atoms) * 2:  # 2 spin states per atom
                spec_mlk = {}
                for i_spec in species:
                    if i_spec not in spec_mlk:
                        spec_mlk[i_spec] = 0
                    i_spec_id = species_id[i_spec]
                    for j in i_spec_id:
                        spec_mlk[i_spec] += currentMlk[j]  # adding up all contributions of each atom type
                maxm = 0
                max_spec = 'Pb'
                for spec in species:
                    if spec_mlk[spec] > maxm:
                        maxm = spec_mlk[spec]
                        max_spec = spec
                plot_dot = False
                if currentState % 2 == 0:
                    if currentK % mod == 0:
                        plot_dot = True
                else:
                    if (currentK + int(mod / 2)) % mod == 0:
                        plot_dot = True
                cur_x = xvals[file_id][currentK - 1]
                cur_E = energys[file_id][currentK - 1][currentState - 1]
                # if 0.5 < cur_E < 2.7 and abs(cur_x - (0.5 * round(cur_x / 0.5))) < 0.05:
                #     print(str(cur_x) + ", " + str(cur_E))
                if plot_dot:  # this if statement exists just so "dots" are not too close together in output
                    plt.plot(cur_x, cur_E, 'o', color=species_color[max_spec], markersize=markersizeunit,
                             markeredgecolor=species_color[max_spec])
            i += 1
        print(str(time.time() - curTime) + " seconds")

    x_pts = [0]
    for kpoint_x in band_len_tot[1:]:
        x_pts.append(kpoint_x)
        plt.axvline(kpoint_x, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
        plt.axhline(0, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
    x_pts.append(tot_len)
    # plt.yticks(range(ymin, ymax + 1), [])
    if labels == 0:
        print("No labels found")
    else:
        plt.xticks(ticks=x_pts, labels=labels)
    plt.ylabel('Energy (eV)')
    plt.xlabel("Wave vector, $k$")
    plt.title(title, weight='bold')
    plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
    plt.tight_layout()
    if save and not filename == 0:
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()
    plt.clf()


def mulliken_plot(settings_file, alt_geo=False, debug=False, quiet=False, save=False):
    ############################################
    # Setup                                    #
    ############################################
    input_settings = read_band_settings(settings_file, debug)
    filepath = input_settings[0]
    filename = input_settings[1]
    energyshift = input_settings[2]
    ymin = input_settings[3]
    ymax = input_settings[4]
    substate = input_settings[5]
    color_dict = input_settings[6]
    labels = input_settings[7]
    title = input_settings[8]
    eq = input_settings[9]
    bands = input_settings[10]
    mod = int(input_settings[11])

    flags = 1  # normal run state, print out timing
    if debug:
        flags = 0  # full debug state, print all info
    elif quiet:
        flags = 2  # quiet run state, print out nothing

    ############################################
    # geometry.in                              #
    ############################################
    if not alt_geo:
        rlatvec, atoms, species, species_id = bg.read_geo_for_bands(filepath, color_dict, flags)
    else:
        rlatvec, atoms, species, species_id = bg.read_geo_for_bands(alt_geo, color_dict, flags)
    ############################################
    # control.in                               #
    ############################################
    xvals, band_len_tot, tot_len = bc.read_control_for_bands(filepath, rlatvec, bands, eq, debug)

    ############################################
    # pre-processing                           #
    ############################################
    species_color = color_dict
    if not species_color:
        python_colors = ['r', 'y', 'b', 'c', 'm', 'g', 'pink', 'purple', 'black']
        for s in range(len(species)):
            species_color[species[s]] = python_colors[s]

    sub_st = False
    if not substate == 0:
        sub_st = True  # 0, s; 1, p; 2, d; ...

    plt.clf()
    linewidth = 1
    black_bands = 1
    markersizeunit = 2
    line_mlk_start_id = 5  # column of contributions in mulliken output
    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    matplotlib.rc('text', usetex='true')
    matplotlib.rcParams['axes.linewidth'] = 2
    matplotlib.rcParams.update({'font.size': 20})

    all_files = [join(filepath, f) for f in listdir(filepath) if isfile(join(filepath, f))]
    band_files = [f for f in all_files if 'band' in f.split("/")[-1] and f.endswith('.out') and 'bandmlk' not in f]
    band_files.sort()
    band_files = [band_files[x] for x in bands]
    band_mlkfiles = [f for f in all_files if 'bandmlk' in f and f.endswith('.out')]
    band_mlkfiles.sort()
    mlk = False
    if len(band_mlkfiles):
        mlk = True
        band_mlkfiles = [band_mlkfiles[x] for x in bands]

    if not flags == 2:
        print("color key: ", species_color)
        if flags == 0:
            print(band_files)
            if mlk:
                print(band_mlkfiles)

    eshift = False
    for shift in energyshift:
        if shift != 0:
            eshift = True

    ############################################
    # band10**.out							   #
    ############################################
    max_occ = -100
    if not flags == 2:
        print("Analyzing band10**.out files...")
    energys = []
    for file_id in range(len(band_files)):
        temp_en = []
        file = band_files[file_id]
        curTime = time.time()
        if not flags == 2:
            print("Processing " + file, end='\t\t')
        f = open(file)
        for line in f:
            words = line.strip().split()
            energy = []
            occ_ene = words[4:]
            for i in range(0, len(occ_ene), 2):
                occupation = occ_ene[i]
                cur_energy = float(occ_ene[i + 1])
                if occupation == "1.00000":
                    if cur_energy > max_occ:
                        max_occ = cur_energy
                energy.append(cur_energy - energyshift[file_id])
            temp_en.append(energy)
        # len(energy) = number of states; len(temp_en) = number of k points
        # band_points.append(len(k_index))
        energys.append(temp_en)
        if not flags == 2:
            print(str(time.time() - curTime) + " seconds")
    if not eshift:
        energys = [[[val - max_occ for val in kpt] for kpt in band] for band in energys]
        if not flags == 2:
            if max_occ > 0:
                print(f'Max occupied state found at %s eV' % max_occ)
                return
            print(f'Max occupied state found at %s eV' % max_occ)

    if not flags == 2:
        print("Plotting bands", end='\t\t')
    curTime = time.time()
    for count, path in enumerate(energys):
        for i in range(len(path[0])):
            band = []
            for j in range(len(path)):
                band.append(path[j][i])
            if ymin - 1 < band[0] < ymax + 1:
                plt.plot(xvals[count], band, color='k', lw=black_bands)
    if not flags == 2:
        print(str(time.time() - curTime) + " seconds")

    ############################################
    # output non-mlk						   #
    ############################################
    if not mlk:
        x_pts = []
        x_pts.append(0)
        for kpoint_x in band_len_tot[1:]:
            x_pts.append(kpoint_x)
            plt.axvline(kpoint_x, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
            plt.axhline(0, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])

        x_pts.append(tot_len)
        if labels == 0:
            print("No labels found")
        else:
            plt.xticks(ticks=x_pts, labels=labels)
        plt.ylabel('Energy (eV)')
        plt.xlabel("Wave vector, $k$")
        plt.title(title, weight='bold')
        plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
        plt.tight_layout()
        if save and not filename == 0:
            plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
        else:
            plt.show()
        plt.clf()
        return 0

    ############################################
    # mlk									   #
    ############################################
    if not flags == 2:
        print("Analyzing bandmlk100*.out files...")

    # for file_id in range(len(band_mlkfiles)):
    for file_id, mlkfile in enumerate(band_mlkfiles):
        mlkfile = band_mlkfiles[file_id]
        curTime = time.time()
        print("Processing " + mlkfile, end='\t\t')
        f = open(mlkfile)
        currentMlk = []
        i = 0
        for line in f:
            if line.strip() == '':
                continue
            elif line.strip().startswith("k point number:"):
                currentK = int(line.strip().split()[3][:-1])
                kpt = line.strip()
                currentMlk = []
            elif line.strip().startswith("State") and len(line.strip().split()) == 2:
                words = line.strip().split()
                currentState = int(words[1])
                currentMlk = []
            elif line.strip()[0].isdigit():
                words = line.strip().split()
                if sub_st:
                    currentMlk.append(float(words[line_mlk_start_id + 1 + substate]))
                else:
                    only_pos = False
                    if only_pos:
                        temp = 0
                        for x in words[line_mlk_start_id + 1:]:
                            if float(x) > 0:
                                temp += float(x)
                        currentMlk.append(temp)
                    else:
                        currentMlk.append(float(words[line_mlk_start_id]))
            current_energy = float(words[1])
            if len(currentMlk) == len(atoms) * 2:  # 2 spin states per atom
                spec_mlk = {}
                for i_spec in species:
                    if i_spec not in spec_mlk:
                        spec_mlk[i_spec] = 0
                    i_spec_id = species_id[i_spec]
                    for j in i_spec_id:
                        spec_mlk[i_spec] += currentMlk[j]  # adding up all contributions of each atom type
                maxm = 0
                max_spec = 'Pb'
                for spec in species:
                    if spec_mlk[spec] > maxm:
                        maxm = spec_mlk[spec]
                        max_spec = spec
                plot_dot = False
                if currentState % 2 == 0:
                    if currentK % mod == 0:
                        plot_dot = True
                else:
                    if (currentK + int(mod / 2)) % mod == 0:
                        plot_dot = True
                cur_x = xvals[file_id][currentK - 1]
                cur_E = energys[file_id][currentK - 1][currentState - 1]
                # if 0.5 < cur_E < 2.7 and abs(cur_x - (0.5 * round(cur_x / 0.5))) < 0.05:
                #     print(str(cur_x) + ", " + str(cur_E))
                if plot_dot:  # this if statement exists just so "dots" are not too close together in output
                    if ymin - 1 < cur_E < ymax + 1:
                        plt.plot(cur_x, cur_E, 'o', color=species_color[max_spec], markersize=markersizeunit,
                                 markeredgecolor=species_color[max_spec])
            i += 1
        print(str(time.time() - curTime) + " seconds")

    x_pts = [0]
    for kpoint_x in band_len_tot[1:]:
        x_pts.append(kpoint_x)
        plt.axvline(kpoint_x, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
        plt.axhline(0, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
    x_pts.append(tot_len)
    if labels == 0:
        print("No labels found")
    else:
        plt.xticks(ticks=x_pts, labels=labels)
    plt.ylabel('Energy (eV)')
    plt.xlabel("Wave vector, $k$")
    plt.title(title, weight='bold')
    plt.axis([0, xvals[-1][len(xvals[-1]) - 1], ymin, ymax])
    # plt.axis([0, xvals[-1][0], ymin, ymax])
    plt.tight_layout()
    if save and not filename == 0:
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()
    plt.clf()


def read_band_settings(settings_file, debug=False):
    with open(settings_file, "r") as f:
        settings = f.readlines()
    substate = 0
    eq = False
    mod = 16
    energyshift = [0]
    filepath = "/".join(settings_file.split("/")[:-1])
    figure_name = filepath + "/" + settings_file.split("/")[-1] + ".png"
    for ln in settings:
        temp = ln.replace("\t", " ").split()
        if len(temp) == 0:
            continue
        if temp[0] == "filepath":
            filepath = temp[1]
        elif temp[0] == "figure_name":
            figure_name = temp[1]
        elif temp[0] == "energyshift":
            energyshift = [float(x) for x in temp[1:]]
        elif temp[0] == "bounds":
            ymin = float(temp[1])
            ymax = float(temp[2])
        elif temp[0] == "substate":
            substate = float(temp[1])
        elif temp[0] == "color_dict":
            color_dict = {}
            for x in range(int((len(temp) - 1) / 2)):
                color_dict.update({temp[2 * x + 1]: temp[2 * x + 2]})
        elif temp[0] == "labels":
            labels = temp[1:]
        elif temp[0] == "title":
            title = " ".join(temp[1:])
        elif temp[0] == "equal":
            eq = True
        elif temp[0] == "bands":
            bands = [int(x) for x in temp[1:]]
        elif temp[0] == "modulo":
            mod = float(temp[1])
        else:
            print(f'keyword not recognized: %s' % ln)
    if len(energyshift) < len(bands):
        energyshift = [energyshift[0] for x in bands]
    return [filepath, figure_name, energyshift, ymin, ymax, substate, color_dict, labels, title, eq, bands, mod]


def new_correlation_plot():
    file = "../../FHI-aims/Double_Perovskites/Figures/correlations/new_IP.txt"
    with open(file, "r") as f:
        IP_raw = []
        for line in f.readlines():
            if "#" not in line:
                IP_raw.append(line.split())
    B_avg = [x[0] for x in IP_raw]
    DB = [x[1] for x in IP_raw]
    DE = [x[2] for x in IP_raw]
    value_dict = {}
    for count, value in enumerate(B_avg):
        key = f'%s_%s' % (value, DB[count])
        if key not in value_dict:
            value_dict.update({key: [float(DE[count])]})
        else:
            arr = value_dict.pop(key)
            value_dict.update({key: arr + [float(DE[count])]})

    reduced_dict = {}
    for key in value_dict.keys():
        reduced_dict.update({key: sum(value_dict[key]) / len(value_dict[key])})

    set_dict = {}
    for key in reduced_dict.keys():
        temp = key.split("_")
        if temp[0] not in set_dict:
            set_dict.update({temp[0]: [[int(temp[1])], [reduced_dict[key]]]})
        else:
            arr = set_dict.pop(temp[0])
            new_arr = [arr[0] + [int(temp[1])], arr[1] + [reduced_dict[key]]]
            set_dict.update({temp[0]: new_arr})

    matplotlib.rc('text', usetex='true')
    color_set = ['rosybrown', 'firebrick', 'red', 'lightsalmon', 'olive', 'green', 'mediumseagreen', 'lime',
                 'teal', 'darkturquoise', 'cyan', 'royalblue', 'rebeccapurple', 'blueviolet', 'fuchsia', 'orchid']
    i = 0
    print(set_dict)
    for key in set_dict:
        current = set_dict[key]
        if len(current[0]) == 1:
            plt.plot(current[0], current[1], 'o', color='k', label="_")
        else:
            plt.plot(current[0], current[1], 'o', color=color_set[i], label=key)
            a, b = np.polyfit(current[0], current[1], 1)
            plt.plot(current[0], [a * x + b for x in current[0]], color=color_set[i], label="_")
            i += 1
    plt.ylabel('$\Delta E\pm$ (eV)')
    plt.xlabel(r'$|\Delta\beta|$')
    plt.legend()
    plt.show()


def correlation_plot():
    folder = "../../FHI-aims/Double_Perovskites/Figures/correlations/"

    with open(folder + "IP.txt", "r") as f:
        IP_raw = [x.split() for x in f.readlines()]
    DB_IP = [float(x[0]) for x in IP_raw]
    Maurer_IP = [float(x[5]) for x in IP_raw]
    Sig_sq_IP = [float(x[8]) for x in IP_raw]
    Delt_d_IP = [float(x[11]) for x in IP_raw]
    SS_1_IP = [float(x[12]) for x in IP_raw]
    SS_2_IP = [float(x[13]) for x in IP_raw]
    SS_3_IP = [float(x[14]) for x in IP_raw]
    SS_4_IP = [float(x[15]) for x in IP_raw]
    SS_5_IP = [float(x[16]) for x in IP_raw]
    SS_6_IP = [float(x[17]) for x in IP_raw]

    with open(folder + "OOP.txt", "r") as f:
        OOP_raw = [x.split() for x in f.readlines()]
    DB_OOP = [float(x[0]) for x in OOP_raw]
    Maurer_OOP = [float(x[5]) for x in OOP_raw]
    Sig_sq_OOP = [float(x[8]) for x in OOP_raw]
    Delt_d_OOP = [float(x[11]) for x in OOP_raw]
    SS_1_OOP = [float(x[12]) for x in OOP_raw]
    SS_2_OOP = [float(x[13]) for x in OOP_raw]
    SS_3_OOP = [float(x[14]) for x in OOP_raw]
    SS_4_OOP = [float(x[15]) for x in OOP_raw]
    SS_5_OOP = [float(x[16]) for x in OOP_raw]
    SS_6_OOP = [float(x[17]) for x in OOP_raw]

    with open(folder + "Maurer.txt", "r") as f:
        M_raw = [x.split() for x in f.readlines()]
    end = 5
    DB_M = [float(x[1]) for x in M_raw[:end]]
    Maurer_M = [float(x[0]) for x in M_raw[:end]]
    Sig_sq_M = [float(x[4]) for x in M_raw[:end]]
    Delt_d_M = [float(x[7]) for x in M_raw[:end]]
    SS_1_M = [float(x[8]) for x in M_raw[:end]]
    SS_2_M = [float(x[9]) for x in M_raw[:end]]
    SS_3_M = [float(x[10]) for x in M_raw[:end]]
    SS_4_M = [float(x[11]) for x in M_raw[:end]]
    SS_5_M = [float(x[12]) for x in M_raw[:end]]
    SS_6_M = [float(x[13]) for x in M_raw[:end]]

    with open(folder + "Random.txt", "r") as f:
        Random_raw = [x.split() for x in f.readlines()]
    DB_R = [float(x[0]) for x in Random_raw]
    Maurer_R = [float(x[5]) for x in Random_raw]
    Sig_sq_R = [float(x[8]) for x in Random_raw]
    Delt_d_R = [float(x[11]) for x in Random_raw]
    SS_1_R = [float(x[12]) for x in Random_raw]
    SS_2_R = [float(x[13]) for x in Random_raw]
    SS_3_R = [float(x[14]) for x in Random_raw]
    SS_4_R = [float(x[15]) for x in Random_raw]
    SS_5_R = [float(x[16]) for x in Random_raw]
    SS_6_R = [float(x[17]) for x in Random_raw]

    matplotlib.rc('text', usetex='true')

    # fig = plt.figure(figsize=(15, 8))
    # gs = fig.add_gridspec(4, 4, hspace=0.2, wspace=0.1)
    # axs = gs.subplots(sharey='row')

    msize = 4
    alph = 0.7
    # these seem to be cases of positive delta beta, negative delta beta, all, and then all while split
    if False:
        start1 = 0
        end1 = 5
        start2 = 9
        end2 = 14
        start3 = 18
        end3 = 23
        start4 = 27
        end4 = 32
    elif False:
        start1 = 4
        end1 = 9
        start2 = 13
        end2 = 18
        start3 = 22
        end3 = 27
        start4 = 31
        end4 = 36
    elif False:
        start1 = 0
        end1 = 9
        start2 = 9
        end2 = 18
        start3 = 18
        end3 = 27
        start4 = 27
        end4 = 36
    else:
        start1 = 0
        end1 = 5
        start2 = 9
        end2 = 14
        start3 = 18
        end3 = 23
        start4 = 27
        end4 = 32

        # plt.plot(DB_IP[start1:end1], SS_2_IP[start1:end1], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start2:end2], SS_2_IP[start2:end2], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start3:end3], SS_2_IP[start3:end3], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start4:end4], SS_2_IP[start4:end4], 'o', color='r', markersize=msize, alpha=alph)
        # temp_domain = DB_IP[start1:end1] + DB_IP[start2:end2] + DB_IP[start3:end3] + DB_IP[start4:end4]
        # temp_range = SS_2_IP[start1:end1] + SS_2_IP[start2:end2] + SS_2_IP[start3:end3] + SS_2_IP[start4:end4]

        SS = [max([SS_1_IP[i], SS_2_IP[i], SS_3_IP[i], SS_4_IP[i], SS_5_IP[i], SS_6_IP[i]]) for i in
              range(len(SS_1_IP))]

        start = start4
        end = end4

        plt.plot(DB_IP[start:end], SS[start:end], 'o', color='r', markersize=msize, alpha=alph, label="_")
        temp_domain = DB_IP[start:end]
        temp_range = SS[start:end]
        a, b = np.polyfit(temp_domain, temp_range, 1)
        plt.plot(temp_domain, [a * x + b for x in temp_domain], color='r', label="negative")

        start1 = 4
        end1 = 9
        start2 = 13
        end2 = 18
        start3 = 22
        end3 = 27
        start4 = 31
        end4 = 36

        start = start4
        end = end4

        # plt.plot(DB_IP[start1:end1], SS_2_IP[start1:end1], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start2:end2], SS_2_IP[start2:end2], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start3:end3], SS_2_IP[start3:end3], 'o', color='r', markersize=msize, alpha=alph)
        # plt.plot(DB_IP[start4:end4], SS_2_IP[start4:end4], 'o', color='r', markersize=msize, alpha=alph)
        # temp_domain = DB_IP[start1:end1] + DB_IP[start2:end2] + DB_IP[start3:end3] + DB_IP[start4:end4]
        # temp_range = SS_2_IP[start1:end1] + SS_2_IP[start2:end2] + SS_2_IP[start3:end3] + SS_2_IP[start4:end4]

        plt.plot(DB_IP[start:end], SS[start:end], 'o', color='b', markersize=msize, alpha=alph, label="_")
        temp_domain = DB_IP[start:end]
        temp_range = SS[start:end]
        a, b = np.polyfit(temp_domain, temp_range, 1)
        plt.plot(temp_domain, [a * x + b for x in temp_domain], color='b', label="positive")
        plt.ylabel('$\Delta E\pm$ (eV)')
        plt.xlabel(r'$|\Delta\beta|$')
        plt.title(r'Base Structure: $150^{\circ}$')
        plt.legend()
        plt.ylim([0, 0.20])
        plt.show()

    # IP graphs
    color_set = ['rosybrown', 'firebrick', 'red', 'lightsalmon']
    xvals = [Delt_d_IP, Sig_sq_IP, DB_IP, Maurer_IP]
    SS = SS_2_IP

    for i in range(4):
        axs[0, i].plot(xvals[i][start1:end1], SS[start1:end1], 'o', color=color_set[i], markersize=msize, alpha=alph)
        axs[0, i].plot(xvals[i][start2:end2], SS[start2:end2], 's', color=color_set[i], markersize=msize, alpha=alph)
        axs[0, i].plot(xvals[i][start3:end3], SS[start3:end3], 'P', color=color_set[i], markersize=msize, alpha=alph)
        axs[0, i].plot(xvals[i][start4:end4], SS[start4:end4], '^', color=color_set[i], markersize=msize, alpha=alph)
        temp_domain = xvals[i][start1:end1] + xvals[i][start2:end2] + xvals[i][start3:end3] + xvals[i][start4:end4]
        temp_range = SS[start1:end1] + SS[start2:end2] + SS[start3:end3] + SS[start4:end4]
        a, b = np.polyfit(temp_domain, temp_range, 1)
        axs[0, i].plot(temp_domain, [a * x + b for x in temp_domain], color=color_set[i])
        r = np.corrcoef(temp_domain, temp_range)[0, 1]
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        axs[0, i].text(min(temp_domain), max(temp_range) * 0.9, f'r=%0.4f' % r, bbox=bbox)

    # OOP graphs
    color_set = ['olive', 'green', 'mediumseagreen', 'lime']
    xvals = [Delt_d_OOP, Sig_sq_OOP, DB_OOP, Maurer_OOP]
    SS = SS_2_OOP

    for i in range(4):
        axs[1, i].plot(xvals[i][start1:end1], SS[start1:end1], 'o', color=color_set[i], markersize=msize, alpha=alph)
        axs[1, i].plot(xvals[i][start2:end2], SS[start2:end2], 's', color=color_set[i], markersize=msize, alpha=alph)
        axs[1, i].plot(xvals[i][start3:end3], SS[start3:end3], 'P', color=color_set[i], markersize=msize, alpha=alph)
        axs[1, i].plot(xvals[i][start4:end4], SS[start4:end4], '^', color=color_set[i], markersize=msize, alpha=alph)
        temp_domain = xvals[i][start1:end1] + xvals[i][start2:end2] + xvals[i][start3:end3] + xvals[i][start4:end4]
        temp_range = SS[start1:end1] + SS[start2:end2] + SS[start3:end3] + SS[start4:end4]
        a, b = np.polyfit(temp_domain, temp_range, 1)
        axs[1, i].plot(temp_domain, [a * x + b for x in temp_domain], color=color_set[i])
        r = np.corrcoef(temp_domain, temp_range)[0, 1]
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        axs[1, i].text(min(temp_domain), max(temp_range) * 0.9, f'r=%0.4f' % r, bbox=bbox)
        # ci = 1.96 * np.std(temp_range)/np.sqrt(len(temp_domain))
        # y = [a * x + b for x in temp_domain]
        # axs[1, i].fill_between(temp_domain, (y - ci), (y + ci), color=color_set[i], alpha=0.1)

    # Maurer graphs
    color_set = ['teal', 'darkturquoise', 'cyan', 'royalblue']
    xvals = [Delt_d_M, Sig_sq_M, DB_M, Maurer_M]
    SS = SS_2_M

    for i in range(4):
        axs[2, i].plot(xvals[i], SS, 'o', color=color_set[i], markersize=msize, alpha=alph)
        temp_domain = xvals[i]
        temp_range = SS
        a, b = np.polyfit(temp_domain, temp_range, 1)
        axs[2, i].plot(temp_domain, [a * x + b for x in temp_domain], color=color_set[i])
        r = np.corrcoef(temp_domain, temp_range)[0, 1]
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        axs[2, i].text(min(temp_domain), max(temp_range) * 0.9, f'r=%0.4f' % r, bbox=bbox)

    # Random graphs
    color_set = ['rebeccapurple', 'blueviolet', 'fuchsia', 'orchid']
    xvals = [Delt_d_R, Sig_sq_R, DB_R, Maurer_R]
    SS = SS_2_R

    for i in range(4):
        axs[3, i].plot(xvals[i], SS, 'o', color=color_set[i], markersize=msize, alpha=alph)
        temp_domain = xvals[i]
        temp_range = SS
        a, b = np.polyfit(temp_domain, temp_range, 1)
        axs[3, i].plot(temp_domain, [a * x + b for x in temp_domain], color=color_set[i])
        r = np.corrcoef(temp_domain, temp_range)[0, 1]
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        axs[3, i].text(min(temp_domain), max(temp_range) * 0.9, f'r=%0.4f' % r, bbox=bbox)

    # labels
    axs[3, 0].set(xlabel='$\Delta d$')
    axs[3, 1].set(xlabel='$\sigma^2$')
    axs[3, 2].set(xlabel='$\Delta\\beta$')
    axs[3, 3].set(xlabel='$d_{diag}$')

    axs[0, 0].set(ylabel='$\Delta E\pm$ (eV)')
    axs[1, 0].set(ylabel='$\Delta E\pm$ (eV)')
    axs[2, 0].set(ylabel='$\Delta E\pm$ (eV)')
    axs[3, 0].set(ylabel='$\Delta E\pm$ (eV)')

    filename = "../../FHI-aims/Double_Perovskites/Figures/images/fig4_all.png"
    # plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    plt.show()


# noinspection PyTypeChecker
def plot_3d_solid_with_path_and_names(geo_file, corners, adjacency, pathway, names, setup=False, save=False,
                                      filename=None):
    # use: https://lampz.tugraz.at/~hadley/ss1/bzones/drawing_BZ.php
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if setup:
        plt.show()
    # Plot pathway
    color_dict = ['r', 'orange', 'y', 'g', 'b']
    color_dict = color_dict * 50
    reciprocal = bg.reciprocal_vectors(geo_file)
    pathway_recip = []
    for i in range(len(pathway)):
        point = pathway[i][0:3]
        temp_x1 = reciprocal[0][0] * point[0] + reciprocal[1][0] * point[1] + reciprocal[2][0] * point[2]
        temp_y1 = reciprocal[0][1] * point[0] + reciprocal[1][1] * point[1] + reciprocal[2][1] * point[2]
        temp_z1 = reciprocal[0][2] * point[0] + reciprocal[1][2] * point[1] + reciprocal[2][2] * point[2]
        point = pathway[i][3:]
        temp_x2 = reciprocal[0][0] * point[0] + reciprocal[1][0] * point[1] + reciprocal[2][0] * point[2]
        temp_y2 = reciprocal[0][1] * point[0] + reciprocal[1][1] * point[1] + reciprocal[2][1] * point[2]
        temp_z2 = reciprocal[0][2] * point[0] + reciprocal[1][2] * point[1] + reciprocal[2][2] * point[2]
        ax.plot([temp_x1, temp_x2],
                [temp_y1, temp_y2],
                [temp_z1, temp_z2], color_dict[i], linewidth=2, marker='o', markersize=7)
        pathway_recip.append([temp_x1, temp_y1, temp_z1, temp_x2, temp_y2, temp_z2])

    for count, corner in enumerate(corners):
        ax.scatter(corner[0], corner[1], corner[2], c='red', marker='o', alpha=0.3)
        # Text for determining adjacency matrix
        if setup:
            ax.text(corner[0], corner[1], corner[2] + 0.05, str(count), color='r', fontsize=14)
        for ind in adjacency[count]:
            ax.plot([corners[ind][0], corner[0]],
                    [corners[ind][1], corner[1]],
                    [corners[ind][2], corner[2]], 'k')

    # Plot names near points
    x_offset = -0.03
    y_offset = 0.02
    z_offset = 0.02
    recip_labels = ["x*", "y*", "z*"]
    for i, recip in enumerate(reciprocal):
        ax.arrow3D(0, 0, 0, 0.75 * recip[0], 0.75 * recip[1], 0.75 * recip[2],
                   mutation_scale=20,
                   arrowstyle="-|>",
                   linestyle='dashed',
                   ec='k',
                   alpha=0.5)
        if i > 0:
            ax.text(0.75 * recip[0] + x_offset, 0.75 * recip[1] + y_offset, 0.75 * recip[2] + z_offset, recip_labels[i],
                    color='k', fontsize=14)
        else:
            ax.text(0.75 * recip[0] - 0.3 * x_offset, 0.75 * recip[1] + y_offset, 0.75 * recip[2] + z_offset,
                    recip_labels[i],
                    color='k', fontsize=14)

    ax.set_axis_off()
    ax.set_proj_type('ortho')
    ax.view_init(elev=65, azim=70, roll=160)
    # plt.show()

    outline_width = 1
    for i, point in enumerate(pathway_recip):
        if sum([abs(x) for x in point[0:3]]) != 0:
            ax.text(point[0] + x_offset, point[1] + y_offset, point[2] + z_offset,
                    names[i][0], color=color_dict[i], fontsize=14,
                    path_effects=[pe.withStroke(linewidth=outline_width, foreground="black")])
        if sum([abs(x) for x in point[3:]]) != 0:
            ax.text(point[3] + x_offset, point[4] + y_offset, point[5] + z_offset,
                    names[i][1], color=color_dict[i], fontsize=14,
                    path_effects=[pe.withStroke(linewidth=outline_width, foreground="black")])

    ax.text(x_offset, y_offset, z_offset, 'Γ', color='k', fontsize=14)
    ax.scatter(0, 0, 0, c='k', marker='o', s=70)

    if save:
        if filename is None:
            filename = "/".join(geo_file.split("/")[:-1]) + "/BZ.png"
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()
    plt.clf()


# histograms for NMR calculations
def NMR_histogram(folders):
    datas = [bao.all_shieldings(folder) for folder in folders]
    # TMS reference is 31.1846 ppm
    # TMS reference new is 31.03126 ppm
    datas = [[x - 31.03126 for x in data] for data in datas]
    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.set_figheight(6)
    fig.subplots_adjust(hspace=0)

    axs[0].hist(datas[0], bins=[-9 + 0.15 * x for x in range(100)], density=True, rwidth=0.8)
    axs[1].hist(datas[1], bins=[-9 + 0.15 * x for x in range(100)], density=True, rwidth=0.8)
    axs[2].hist(datas[2], bins=[-9 + 0.15 * x for x in range(100)], density=True, rwidth=0.8)

    plt.xlim([-9, 0.5])
    plt.xticks(range(-9, 1))
    plt.xlabel("Shieldings (ppm)")
    axs[0].set_ylabel("TA")
    axs[1].set_ylabel("TA(I)3")
    axs[2].set_ylabel("TA(I3)3")

    axs[0].set_yticks([])
    axs[1].set_yticks([])
    axs[2].set_yticks([])

    axs[0].grid(axis='x')
    axs[1].grid(axis='x')
    axs[2].grid(axis='x')

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['bottom'].set(linewidth=2)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].spines['left'].set_visible(False)
    axs[1].spines['bottom'].set(linewidth=2)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    axs[2].spines['left'].set_visible(False)
    plt.show()


def plot_chem_book(width=0.001):
    fig, ax = plt.subplots()
    signals = [(3.811, 109), (3.730, 427), (3.652, 469), (3.576, 132), (2.607, 199), (2.599, 168),
               (1.303, 430), (1.286, 23), (1.226, 1000), (1.207, 36), (1.199, 25), (1.146, 376)]
    red = -4.087999375000001
    blue = -2.6221356250000003
    green = -1.4276929166666676

    offset = 0.0725 / 2
    # signals = [(red + 3 * offset, 4), (red + offset, 12), (red - offset, 12), (red - 3 * offset, 4),
    #            (blue, 3), (blue, 3),
    #            (green + 2 * offset, 9), (green, 18), (green - 2 * offset, 9)]
    x_range = np.arange(-5, 0.5, 0.0005)
    red_total = [0 for x in x_range]
    blue_total = [0 for x in x_range]
    green_total = [0 for x in x_range]
    for count, signal in enumerate(signals):
        if count < 4:
            red_total = [red_total[i] + signal[1] * norm.pdf(x, -1 * signal[0], width) for i, x in enumerate(x_range)]
        elif count < 6:
            blue_total = [blue_total[i] + signal[1] * norm.pdf(x, -1 * signal[0], width) for i, x in enumerate(x_range)]
        else:
            green_total = [green_total[i] + signal[1] * norm.pdf(x, -1 * signal[0], width) for i, x in
                           enumerate(x_range)]

    ax.plot(x_range, red_total, color='r')
    ax.fill_between(x_range, red_total, 0, color='r')

    ax.plot(x_range, blue_total, color='b')
    ax.fill_between(x_range, blue_total, 0, color='b')

    ax.plot(x_range, green_total, color='g')
    ax.fill_between(x_range, green_total, 0, color='g')

    ax.plot(x_range, [0 for x in x_range], color='k')
    plt.xticks([x for x in range(-10, 1)], [20, 40, 60, 80, 100, 5, 4, 3, 2, 1, 0])
    plt.ylim([0, None])
    plt.yticks([])
    plt.xlim([-5, 0.5])
    plt.xlabel("ppm")
    for spine in ax.spines.values():
        spine.set_linewidth(3)  # Adjust the thickness here
    # plt.box(False)
    plt.tight_layout()
    plt.show()


def NMR_density(folder, atom_dict=None, color_dict=None, average=False, MD=False, width=0.01):
    fig, ax = plt.subplots()
    if MD:
        shieldings = bao.all_shieldings(folder)
    else:
        shieldings = bao.NMR_shielding_values(folder + "aims.out")

    print(shieldings)

    atoms = [x[0] for x in shieldings]
    shields = [x[1] for x in shieldings]
    x_range = np.arange(min(shields) - 31.03126 - .2, max(shields) - 31.03126 + 0.2, 0.01)
    x_range = np.arange(-5, 0.5, 0.0001)
    # x_range = np.arange(-5, 0.5, 0.01)

    if atom_dict is None:
        atom_dict = [atoms]

    total = [0 for x in x_range]
    if not color_dict:
        color_dict = ['r', 'g', 'b', 'purple']
        color_dict = ['lightgreen', 'b', 'r']
    for count, atom_set in enumerate(atom_dict):
        y_values = [0 for x in x_range]
        centers = []
        for pair in shieldings:
            if pair[0] in atom_set:
                if not average:
                    temp = [norm.pdf(x, pair[1] - 31.1846, width) for x in x_range]
                    total = [total[x] + temp[x] for x in range(len(temp))]
                    y_values = [y_values[x] + temp[x] for x in range(len(temp))]
                else:
                    centers.append(pair[1] - 31.1846)
        if average:
            val = sum(centers) / len(centers)
            # print(val)
            temp = [norm.pdf(x, val, width) for x in x_range]
            total = [total[x] + temp[x] for x in range(len(temp))]
            temp = [len(centers) * x for x in temp]
            y_values = [y_values[x] + temp[x] for x in range(len(temp))]
        # y_values = [val if val > 0.0001 else -1 for val in y_values]
        ax.plot(x_range, y_values, color=color_dict[count])
        ax.fill_between(x_range, y_values, 0, color=color_dict[count])
    ax.plot(x_range, [0 for x in x_range], color='k')
    plt.xticks([x for x in range(-10, 1)], [20, 40, 60, 80, 100, 5, 4, 3, 2, 1, 0])
    plt.xlim([-5, 0.5])
    plt.ylim([0, None])
    plt.yticks([])
    plt.xlabel("ppm")
    for spine in ax.spines.values():
        spine.set_linewidth(3)  # Adjust the thickness here
    plt.tight_layout()
    plt.show()


def NMR_average(folder, atom_dict, color_dict=None, width=0.01, xlim=(-9.5, 0.5), type='H'):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    data = {}
    for direct in all_dir:
        if os.path.isdir(folder + direct):
            shields = folder + direct + "/shieldings.out"
            if not os.path.exists(shields):
                bao.create_shield_out(folder + direct + "/")
            with open(shields, "r") as f:
                lines = f.readlines()
                if len(lines) == 3:
                    continue
                for line in lines:
                    if "#" not in line:
                        temp = line.split()
                        if int(temp[0]) not in data:
                            data.update({int(temp[0]): [float(temp[2])]})
                        else:
                            data[int(temp[0])].append(float(temp[2]))

    for key in data:
        data[key] = sum(data[key]) / len(data[key])
        # data[key] = data[key] - 31.1846
        if type == 'H':
            data[key] = data[key] - 31.03126
        elif type == 'C':
            continue
            data[key] = data[key] - 177.4423
        elif type == 'N':
            # data[key] = data[key] - 257.1237
            data[key] = data[key] + 143.3070
        else:
            print("Unrecognized species passed into 'type' parameter")
            return

    x_range = np.arange(xlim[0], xlim[1], 0.001)
    if not color_dict:
        color_dict = ['r', 'c', 'b', 'g']
        color_dict = ['r', 'g', 'b']

    for count, structure_set in enumerate(atom_dict):
        y_values = [0 for x in x_range]
        for atom_set in structure_set:
            centers = []
            for atom in atom_set:
                centers.append(data[atom])
            center = sum(centers) / len(centers)
            print(atom_set)
            print(centers)
            print(center)
            temp = [norm.pdf(x, center, width) for x in x_range]
            factor = math.sqrt(len(centers))
            temp = [factor * x for x in temp]
            y_values = [y_values[x] + temp[x] for x in range(len(temp))]
        plt.fill_between(x_range, y_values, color=color_dict[count])
    plt.plot(x_range, [0 for x in x_range], color='k')
    ticks = [x for x in range(int(xlim[0] - 1), int(xlim[1]) + 1)]
    # print(ticks)
    # while len(ticks) > 20:
    #     ticks = [ticks[2 * x] for x in range(int(len(ticks) / 2))]
    plt.xticks(ticks)
    plt.yticks([])
    plt.xlim(xlim)
    plt.xlabel("ppm")
    plt.show()


def visualize_spin_texture_directions(file):
    points = []
    with open(file, "r") as f:
        for line in f:
            if "band_mulliken" in line:
                temp = line.split()
                points.append([float(x) for x in temp[2:8]])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print(points)
    for point in points:
        ax.plot([point[0], point[3]],
                [point[1], point[4]],
                [point[2], point[5]], 'k')
    plt.show()


def plot_beta_avg_corr(data_file):
    matplotlib.rc('text', usetex='true')
    data = []
    with open(data_file, "r") as f:
        for line in f.readlines():
            if "#" not in line:
                data.append(line.split())
    delta_beta = [float(x[0]) for x in data]
    data_180_Y = [float(x[16]) for x in data]
    data_170_Y = [float(x[11]) for x in data]
    data_160_Y = [float(x[6]) for x in data]
    data_150_Y = [float(x[1]) for x in data]
    current_data = [data_180_Y, data_170_Y, data_160_Y, data_150_Y]
    colors = ['r', 'g', 'b', 'y']
    labels = ["180", "170", "160", "150"]
    for count, series in enumerate(current_data):
        plt.plot(delta_beta, series, color=colors[count], label=labels[count])
        plt.scatter(delta_beta, series, color=colors[count], label="_")
    plt.legend()
    plt.xlabel(r'$$\Delta\beta (^{\circ})$$')
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $Y\rightarrow\Gamma$')
    plt.ylim([0, 0.4])
    plt.show()

    data_180_X = [float(x[17]) for x in data]
    data_170_X = [float(x[12]) for x in data]
    data_160_X = [float(x[7]) for x in data]
    data_150_X = [float(x[2]) for x in data]
    current_data = [data_180_X, data_170_X, data_160_X, data_150_X]
    colors = ['r', 'g', 'b', 'y']
    labels = ["180", "170", "160", "150"]
    for count, series in enumerate(current_data):
        plt.plot(delta_beta, series, color=colors[count], label=labels[count])
        plt.scatter(delta_beta, series, color=colors[count], label="_")
    plt.legend()
    plt.xlabel(r'$$\Delta\beta (^{\circ})$$')
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $\Gamma\rightarrow X$')
    plt.ylim([0, 0.4])
    plt.show()

    data_180_M = [float(x[18]) for x in data]
    data_170_M = [float(x[13]) for x in data]
    data_160_M = [float(x[8]) for x in data]
    data_150_M = [float(x[3]) for x in data]
    current_data = [data_180_M, data_170_M, data_160_M, data_150_M]
    colors = ['r', 'g', 'b', 'y']
    labels = ["180", "170", "160", "150"]
    for count, series in enumerate(current_data):
        plt.plot(delta_beta, series, color=colors[count], label=labels[count])
        plt.scatter(delta_beta, series, color=colors[count], label="_")
    plt.legend()
    plt.xlabel(r'$$\Delta\beta (^{\circ})$$')
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $M\rightarrow\Gamma$')
    plt.ylim([0, 0.4])
    plt.show()

    data_180_P = [float(x[19]) for x in data]
    data_170_P = [float(x[14]) for x in data]
    data_160_P = [float(x[9]) for x in data]
    data_150_P = [float(x[4]) for x in data]
    current_data = [data_180_P, data_170_P, data_160_P, data_150_P]
    colors = ['r', 'g', 'b', 'y']
    labels = ["180", "170", "160", "150"]
    for count, series in enumerate(current_data):
        plt.plot(delta_beta, series, color=colors[count], label=labels[count])
        plt.scatter(delta_beta, series, color=colors[count], label="_")
    plt.legend()
    plt.xlabel(r'$$\Delta\beta (^{\circ})$$')
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $\Gamma\rightarrow P$')
    plt.ylim([0, 0.4])
    plt.show()


def plot_correlations(folder, mode="dd", band_path=1, title="put in a title"):
    matplotlib.rc('text', usetex='true')
    base_structures = ["180/", "170/", "160/", "150/"]
    colors = ['r', 'g', 'b', 'y']
    labels = ["180", "170", "160", "150"]
    for count, base in enumerate(base_structures):
        print(base)
        cur_x_vals = []
        if mode == "db":
            cur_x_vals = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20]
        cur_y_vals = []
        all_dir = ['00.0', '02.5', '05.0', '07.5', '10.0', '12.5', '15.0', '17.5', '20.0']
        for directory in all_dir:
            print(directory)
            if mode == "dd":
                cur_x_vals.append(Geometry.new_delt_d(folder + base + directory + "/geometry.in"))
            elif mode == "ss":
                cur_x_vals.append(float(Geometry.new_sigma_sq(folder + base + directory + "/geometry.in")))
            cur_y_vals.append(float(bbo.band_info(folder + base + directory + "/", f'band100{band_path}.out',
                                                  band_gap=False, spin_splitting=True, verbose=False)))

        plt.plot(cur_x_vals, cur_y_vals, color=colors[count], label=labels[count])
        plt.scatter(cur_x_vals, cur_y_vals, color=colors[count], label="_")
    plt.legend()
    if mode == "db":
        plt.xlabel(r'$$\Delta\beta (^{\circ})$$')
    elif mode == "dd":
        plt.xlabel(r'$$\Delta d$$')
    elif mode == "ss":
        plt.xlabel(r'$$\sigma^2$$')
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(title)
    plt.ylim([-0.01, 0.40])
    plt.show()


def plot_maurer_structures(folders, path, mode='ddiag'):
    matplotlib.rc('text', usetex='true')
    colors = ['b', 'orange', 'g']
    labels = ["AgBi", "TlBi", "Pb"]
    for count, folder in enumerate(folders):
        maurer_x, all_splits, _, _, xlab = get_maurer_data(folder, mode, labels[count])
        plt.plot(maurer_x[:10], all_splits[path][:10], color=colors[count], label=labels[count])
        plt.scatter(maurer_x[:10], all_splits[path][:10], color=colors[count], label="_")
    plt.xlabel(xlab)
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    if path == 0:
        plt.title(r'Band Path: $Y\rightarrow\Gamma$')
    elif path == 1:
        plt.title(r'Band Path: $X\rightarrow\Gamma$')
    elif path == 2:
        plt.title(r'Band Path: $M\rightarrow\Gamma$')
    elif path == 3:
        plt.title(r'Band Path: $P\rightarrow\Gamma$')
    plt.legend()
    plt.show()


def plot_maurer_structure(folder, x_axis="ddiag", type="AgBi"):
    matplotlib.rc('text', usetex='true')
    maurer_x, all_splits, db_x, db_splits, xlab = get_maurer_data(folder, x_axis, type)
    # band path 1
    plt.xlabel(xlab)
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $Y\rightarrow\Gamma$')
    plt.scatter(maurer_x, all_splits[0], color='r')
    plt.scatter(db_x, db_splits[0], color='k')
    plt.legend(["Maurer data", r'$\Delta\beta$ structure data'])
    plt.show()

    # band path 2
    plt.xlabel(xlab)
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $X\rightarrow\Gamma$')
    plt.scatter(maurer_x, all_splits[1], color='r')
    plt.scatter(db_x, db_splits[1], color='k')
    plt.legend(["Maurer data", r'$\Delta\beta$ structure data'])
    plt.show()

    # band path 3
    plt.xlabel(xlab)
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $M\rightarrow\Gamma$')
    plt.scatter(maurer_x, all_splits[2], color='r')
    plt.scatter(db_x, db_splits[2], color='k')
    plt.legend(["Maurer data", r'$\Delta\beta$ structure data'])
    plt.show()

    # band path 4
    plt.xlabel(xlab)
    plt.ylabel(r'$\Delta E\pm$ (eV)')
    plt.title(r'Band Path: $P\rightarrow\Gamma$')
    plt.scatter(maurer_x, all_splits[3], color='r')
    plt.scatter(db_x, db_splits[3], color='k')
    plt.legend(["Maurer data", r'$\Delta\beta$ structure data'])
    plt.show()


def get_maurer_data(folder, x_axis="ddiag", type="AgBi"):
    write_base = folder
    all_splits = [[], [], [], []]
    dbs = []
    dds = []
    sigmas = []
    for i in range(21):
        band_base = write_base + "s" + str(i) + "/"
        for j in range(1, 5):
            all_splits[j - 1].append(float(bbo.band_info(band_base, f'band100%s.out' % j,
                                                         band_gap=False, spin_splitting=True, verbose=False)))
        angles = [bg.angle_info(band_base + "geometry.in", (1, 7, 4), [[0, 0, 0], [0, 0, 0], [0, 0, 0]])[0],
                  bg.angle_info(band_base + "geometry.in", (4, 10, 1), [[0, 0, 0], [0, 0, 0], [1, 0, 0]])[0],
                  360 - bg.angle_info(band_base + "geometry.in", (4, 8, 1), [[0, 0, 0], [0, 0, 0], [0, 1, 0]])[0],
                  360 - bg.angle_info(band_base + "geometry.in", (4, 9, 1), [[0, 0, 0], [0, 0, 0], [1, 1, 0]])[0]]
        dbs.append(abs(np.average(angles[0:2]) - np.average(angles[2:])))
        dds.append(Geometry.new_delt_d(band_base + "geometry.in"))
        sigmas.append(float(Geometry.new_sigma_sq(band_base + "geometry.in")))

    if x_axis == "db":
        maurer_x = dbs
        db_x = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]
        xlab = r'$$\Delta\beta (^{\circ})$$'
    elif x_axis == "dd":
        maurer_x = dds
        xlab = r'$$\Delta d$$'
        if type == "AgBi":
            db_x = [0.0005514869028169361, 0.0005501529081437721, 0.0005461605279851626,
                    0.0005395385800108273, 0.0005303351139641667, 0.000518617442115841,
                    0.0005044721819242762, 0.00048800531092865903, 0.00046934223390719524]
        elif type == "Pb":
            db_x = [0.00027677247337772547, 0.0002758293590776476, 0.0002730095569457845,
                    0.0002683416952875836, 0.00026187350794989306, 0.00025367186433927324,
                    0.0002438228114710582, 0.000232431628072491, 0.00021962289077003672]
        elif type == "TlBi":
            db_x = [0.001115427028095635, 0.0011136361190813597, 0.0011082729916257015,
                    0.001099366450796213, 0.0010869645260631132, 0.0010711345028494134,
                    0.0010519629667287461, 0.00102955586029827, 0.0010040385527622653]
    elif x_axis == "ddiag":
        maurer_x = [i * math.sqrt(2) * 0.01 for i in range(21)]
        xlab = r'$$\mbox{Displacement Proportion: }(d_{diag}/d_{MH})$$'
        db_x = [0 for i in range(21)]
    elif x_axis == "ss":
        maurer_x = sigmas
        xlab = r'$$\sigma^2$$'
        if type == "AgBi":
            db_x = [0.0, 0.017054138260774294, 0.06821612824867226, 0.1534846954767971, 0.27285771545426807,
                    0.4263322131655902, 0.6139043623387497, 0.8355694844988036, 1.0913220478035366]
        elif type == "Pb":
            db_x = [0.0, 0.01687883704959059, 0.067515093119297, 0.1519080027357202, 0.27005628932233194,
                    0.42195816401746966, 0.6076113240167413, 0.8270129504381312, 1.0801597057064776]
        elif type == "TlBi":
            db_x = [0.0, 0.017001365867382292, 0.0680050909942368, 0.15301005780769691, 0.2720144031951563,
                    0.425015517766429, 0.61201004481777, 0.8329938789952098, 1.0879621646535345]
    if type == "AgBi":
        y_sp = [0.000000000, 0.001978047, 0.007910668, 0.017736152, 0.033331424,
                0.048738543, 0.068330372, 0.090566254, 0.121037926]
        x_sp = [0.000000000, 0.001984659, 0.007954340, 0.017762658, 0.033311418,
                0.048747455, 0.068220740, 0.090601591, 0.121028108]
        m_sp = [0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,
                0.000000000, 0.000000000, 0.000000000, 0.000000000]
        p_sp = [0.000000000, 0.003996148, 0.015753161, 0.035129993, 0.064213747,
                0.091035655, 0.130201712, 0.173199313, 0.218628156]
    elif type == "Pb":
        y_sp = [0.000000000, 0.003651229, 0.014492939, 0.032241702, 0.056427182,
                0.086511194, 0.122028863, 0.162165650, 0.206240499]
        x_sp = [0.000000000, 0.003651221, 0.014462439, 0.032233604, 0.056410237,
                0.086589006, 0.121991426, 0.162074105, 0.206117388]
        m_sp = [0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,
                0.000000000, 0.000000000, 0.000000000, 0.000000000]
        p_sp = [0.000000000, 0.007286301, 0.028664433, 0.063268560, 0.109288709,
                0.165281619, 0.229255956, 0.299524430, 0.374601416]
    elif type == "TlBi":
        y_sp = [0.000000000, 0.005262417, 0.020270133, 0.043109369, 0.071773425,
                0.104614053, 0.140621992, 0.179258340, 0.219193419]
        x_sp = [0.000000000, 0.005292400, 0.020301555, 0.043137159, 0.071756205,
                0.104624037, 0.140624206, 0.179268290, 0.219398904]
        m_sp = [0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,
                0.000000000, 0.000000000, 0.000000000, 0.000000000]
        p_sp = [0.000000000, 0.010452739, 0.039593492, 0.082461978, 0.134546223,
                0.192984300, 0.255323007, 0.318973109, 0.382150808]
    return maurer_x, all_splits, db_x, [y_sp, x_sp, m_sp, p_sp], xlab


# Possible modes:
# db_avg          - delta beta average
# db_max          - delta beta max
# disp_diag_ag    - d_diag displacement of Ag octahedra
# disp_diag_bi    - d_diag displacement of Bi octahedra
# disp_diag_avg   - d_diag average
# disp_diag_max   - d_diag max
# disp_avg        - d_disp average
# disp_max        - d_disp max
# ss_avg          - sigma squared average
# ss_max          - sigma squared max
# dd_avg          - delta d average
# dd_max          - delta d max
def plot_random_structures(folders, title="Random Structures"):
    matplotlib.rc('text', usetex='true')
    mode_dict = {"db_avg": r'$\Delta\beta_{avg}$', "db_max": r'$\Delta\beta_{max}$', "delta_db": r'$\Delta\Delta\beta$',
                 "disp_diag_ag": r'$d_{diag,Ag}$', "disp_diag_bi": r'$d_{diag,Bi}$', "disp_diag_avg": r'$d_{diag,avg}$',
                 "disp_diag_max": r'$d_{diag,max}$', "disp_avg": r'$d_{disp,avg}$', "disp_max": r'$d_{disp,max}$',
                 "ss_avg": r'$\sigma^2_{avg}$', "ss_max": r'$\sigma^2_{max}$',
                 "dd_avg": r'$\Delta d_{avg}$', "dd_max": r'$\Delta d_{max}$',
                 "dl_avg": r'$\Delta L_{avg}$', "dl_max": r'$\Delta L_{max}$', "L2": r'$L^2$', "Ag_len": r'$Ag_{len}$'}
    for count, folder in enumerate(folders):
        all_x, all_splits = get_random_data(folder)
        colors = ['r', 'g', 'b', 'orange']

    # code for 2D intensity maps
    dbs = ["db_avg", "db_max"]
    dbs = ["db_max"]
    dls = ["dl_avg", "dl_max", "L2"]
    # dls = ["Ag_len"]
    average_splits = [max([all_splits[i][j] for i in range(4)]) for j in range(len(all_splits[0]))]
    for count, split in enumerate(average_splits):
        if all_x["db_max"][count] > 13:
            print(count)
            print(all_x["db_max"][count])
            print(split)
    # for path in range(4):
    for db in dbs:
        for dl in dls:
            # plt.scatter(all_x[db], all_x[dl], c=all_splits[path], cmap='jet', edgecolors='k')
            plt.scatter(all_x[db], all_x[dl], c=average_splits, cmap='jet', edgecolors='k')
            plt.colorbar(label='Max Spin Splitting (eV)')
            plt.xlabel(mode_dict[db])
            plt.ylabel(mode_dict[dl])
            # plt.title(f'Path {path + 1}')
            plt.title("Max Spin Splitting")
            plt.show()
            plt.clf()

    for key in all_x.keys():
        x_data = all_x[key]
        maxes = [max([all_splits[j][i] for j in range(len(all_splits))]) for i in range(len(all_splits[0]))]
        slope, intercept, r_value, _, _ = scipy.stats.linregress(x_data, maxes)
        plt.plot([min(x_data), max(x_data)], [x * slope + intercept for x in [min(x_data), max(x_data)]],
                 color='r', label=f'max spin splitting, $r^2$={r_value ** 2}')
        plt.scatter(x_data, maxes, color='r', label=f'_')
        plt.xlabel(mode_dict[key])
        plt.ylabel(r'$\Delta E\pm$ (eV)')
        plt.legend()
        plt.title(title)
        plt.show()
        plt.clf()

        for index, direction in enumerate(all_splits):
            slope, intercept, r_value, _, _ = scipy.stats.linregress(x_data, direction)
            plt.plot([min(x_data), max(x_data)], [x * slope + intercept for x in [min(x_data), max(x_data)]],
                     color=colors[index], label=f'path {index + 1}, $r^2$={r_value ** 2}')
            plt.scatter(x_data, direction, color=colors[index], label="_")
        plt.xlabel(mode_dict[key])
        plt.ylabel(r'$\Delta E\pm$ (eV)')
        plt.legend()
        plt.title(title)
        # filename = folders[0] + key + ".png"
        # plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
        plt.show()
        plt.clf()


def get_random_data(folder):
    write_base = folder
    all_splits = [[], [], [], []]
    dbs = [[], []]
    dds = []
    ddiags = []
    sigmas = []
    dls = []
    l2 = []
    ag_lens = []
    for i in ["05", "10", "15", "20", "25"]:  # for Pb
        # for i in ["05", "10", "15", "20"]: # for AgBi
        # for i in ["25"]:
        for j in range(1, 11):
            band_base = write_base + i + "/" + str(j).rjust(2, "0") + "/"
            for k in range(1, 5):
                all_splits[k - 1].append(float(bbo.band_info(band_base, f'band100%s.out' % k, band_gap=False,
                                                             spin_splitting=True, verbose=False)))

            angles = bg.random_angle_info(band_base + "geometry.in", method="flat")
            dbs[0].append(abs(np.average([angles[0], angles[1]]) - np.average([angles[2], angles[3]])))
            dbs[1].append(abs(np.average([angles[0], angles[2]]) - np.average([angles[1], angles[3]])))
            dds.append(Geometry.new_delt_d(band_base + "geometry.in"))
            sigmas.append(Geometry.new_sigma_sq(band_base + "geometry.in"))
            ddiags.append(Geometry.ddiag_val(band_base + "geometry.in"))
            dls.append(bg.Delta_L(band_base + "geometry.in"))
            l2.append(np.average(bg.Delta_L(band_base + "geometry.in", mode="L2")))
            ag_lens.append(bg.Ag_len(band_base + "geometry.in"))
            # print(f'{all_splits[0][-1]} {all_splits[1][-1]} {all_splits[2][-1]} {all_splits[3][-1]}')
            # mulliken_plot(f'{band_base}settings_final.in', quiet=True)

    all_x = {"db_avg": [np.average([dbs[0][i], dbs[1][i]]) for i in range(len(dbs[0]))],
             "db_max": [max([dbs[0][i], dbs[1][i]]) for i in range(len(dbs[0]))],
             "disp_diag_ag": [max(data_set[4:]) for data_set in ddiags],
             "disp_diag_bi": [max([data_set[2], data_set[3]]) for data_set in ddiags],
             "disp_diag_avg": [np.average(data_set[2:]) for data_set in ddiags],
             "disp_diag_max": [max(data_set[2:]) for data_set in ddiags],
             "disp_avg": [(data_set[0] + data_set[1]) / 2 for data_set in ddiags],
             "disp_max": [max(data_set[0], data_set[1]) for data_set in ddiags],
             "ss_avg": [(sigma[0] + sigma[1]) / 2 for sigma in sigmas],
             "ss_max": [max(sigma[0], sigma[1]) for sigma in sigmas], "dd_avg": [(dd[0] + dd[1]) / 2 for dd in dds],
             "dd_max": [max(dd) for dd in dds],
             "delta_db": [abs(dbs[0][i] - dbs[1][i]) for i in range(len(dbs[0]))],
             "dl_avg": [sum(data_set) / 4 for data_set in dls],
             "dl_max": [max(data_set) for data_set in dls],
             "L2": l2,
             "Ag_len": ag_lens}

    return all_x, all_splits
