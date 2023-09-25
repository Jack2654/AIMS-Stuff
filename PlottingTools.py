import matplotlib
import matplotlib.pyplot as plt
import math
from os import listdir
from os.path import isfile, join
import numpy as np
import BasicGeo as bg
import BasicControl as bc
import time


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

def mulliken_plot(filepath, filename=0, energyshift=0, ymin=-5, ymax=5, substate=0, color_dict={}, labels=0,
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
                energy.append(float(occ_ene[i + 1]) - energyshift)
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
    # plt.show()
    if not filename == 0:
        plt.savefig(filename, dpi=1000, bbox_inches='tight', format="png")
    plt.clf()


def dos_plot(files):
    for file in files:
        lines = []
        with open(file, "r") as f:
            for ln in f:
                if "#" not in ln:
                    lines.append([float(x) for x in ln.split()[:2]])
        energy = []
        dos = []
        for line in lines:
            energy.append(line[0])
            dos.append(line[1])
        print(file)
        plt.plot(energy, dos)
    plt.show()


def MD_plot(file):
    lines = []
    with open(file, "r") as f:
        for ln in f:
            if "#" not in ln:
                lines.append(ln)
    times = []
    temps = []
    for ln in lines:
        temp = ln.split()
        times.append(float(temp[0]))
        temps.append(float(temp[1]))
    plt.plot(times, temps)
    plt.plot((0, 5), (150, 150))
    plt.xlabel("Time (ps)")
    plt.ylabel("Temperature (K)")
    plt.show()
