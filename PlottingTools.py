import matplotlib
import matplotlib.pyplot as plt
import math
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import BasicGeo as bg
import BasicControl as bc
import BasicAimsOut as bao
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
import matplotlib.patheffects as pe


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


def mulliken_plot_old(filepath, filename=0, energyshift=0, ymin=-5, ymax=5, substate=0, color_dict={}, labels=0,
                      title="Default Title", eq=False, debug=False, quiet=False):
    ############################################
    # Setup                                    #
    ############################################
    flags = 1  # normal run state, print out timing
    if debug:
        flags = 0  # full debug state, print all info
    elif quiet:
        flags = 2  # quiet run state, print out nothing

    ############################################
    # geometry.in                              #
    ############################################
    rlatvec, atoms, species, species_id = bg.read_geo_for_bands(filepath, color_dict, flags)

    ############################################
    # control.in                               #
    ############################################
    xvals, band_len_tot, tot_len = bc.read_control_for_bands(filepath, rlatvec, eq, debug)

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
    band_files = [f for f in all_files if 'band' in f and f.endswith('.out') and 'bandmlk' not in f]
    band_files.sort()
    band_mlkfiles = [f for f in all_files if 'bandmlk' in f and f.endswith('.out')]
    band_mlkfiles.sort()
    mlk = False
    if len(band_mlkfiles):
        mlk = True

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
        # x_pts.append(all)
        # plt.yticks(range(ymin, ymax + 1), [])
        plt.xticks(ticks=x_pts, labels=('$X$', '$\Gamma$', '$Y\|L$', '$\Gamma$', '$K$', 'a', 'b'))
        plt.ylabel('Energy (eV)')
        plt.xlabel("Wave vector, $k$")
        # plt.title("3.9\% Doping", weight='bold')
        plt.title("(S-BrMBA)$_2$PbI$_4$", weight='bold')
        plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
        plt.tight_layout()
        plt.show()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
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
                    if currentK % 16 == 0:
                        plot_dot = True
                else:
                    if (currentK + 8) % 16 == 0:
                        plot_dot = True
                cur_x = xvals[file_id][currentK - 1]
                cur_E = energys[file_id][currentK - 1][currentState - 1]
                # if 0.5 < cur_E < 2.7 and abs(cur_x - (0.5 * round(cur_x / 0.5))) < 0.05:
                #     print(str(cur_x) + ", " + str(cur_E))
                if plot_dot:  # this if statement exists just so "dots" are not too close together in output
                    plt.plot(cur_x, cur_E, species_color[max_spec] + 'o', markersize=markersizeunit,
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
    plt.show()
    # if not filename == 0:
    # plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    plt.clf()


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
    colorkey = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
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
            if element == "Pb" or element == "I":
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
        for ln in f:
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


def read_band_settings(settings_file, debug=False):
    with open(settings_file, "r") as f:
        settings = f.readlines()
    substate = 0
    eq = False
    mod = 16
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

    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(4, 4, hspace=0.2, wspace=0.1)
    axs = gs.subplots(sharey='row')

    msize = 4
    alph = 0.7
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
    else:
        start1 = 0
        end1 = 9
        start2 = 9
        end2 = 18
        start3 = 18
        end3 = 27
        start4 = 27
        end4 = 36

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
    plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    # plt.show()


# noinspection PyTypeChecker
def plot_3d_solid_with_path_and_names(geo_file, corners, adjacency, pathway, names, setup=False, save=False, filename=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for count, corner in enumerate(corners):
        ax.scatter(corner[0], corner[1], corner[2], c='red', marker='o', alpha=0.3)
        # Text for determining adjacency matrix
        if setup:
            ax.text(corner[0], corner[1], corner[2] + 0.05, str(count), color='r', fontsize=14)
        for ind in adjacency[count]:
            ax.plot([corners[ind][0], corner[0]],
                    [corners[ind][1], corner[1]],
                    [corners[ind][2], corner[2]], 'k')
    if setup:
        plt.show()
    # Plot pathway
    color_dict = ['r', 'orange', 'y', 'g', 'b']
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

    # Plot names near points
    x_offset = -0.03
    y_offset = 0.02
    z_offset = 0.02
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
            ax.text(0.75 * recip[0] - 0.3 * x_offset, 0.75 * recip[1] + y_offset, 0.75 * recip[2] + z_offset, recip_labels[i],
                    color='k', fontsize=14)

    ax.set_axis_off()
    ax.set_proj_type('ortho')
    ax.view_init(elev=65, azim=70, roll=160)
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
    datas = [[x - 31.1846 for x in data] for data in datas]
    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.set_figheight(6)
    fig.subplots_adjust(hspace=0)

    axs[0].hist(datas[0], bins=[-9 + 0.15*x for x in range(100)], density=True, rwidth=0.8)
    axs[1].hist(datas[1], bins=[-9 + 0.15*x for x in range(100)], density=True, rwidth=0.8)
    axs[2].hist(datas[2], bins=[-9 + 0.15*x for x in range(100)], density=True, rwidth=0.8)

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
