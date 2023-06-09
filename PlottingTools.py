import matplotlib
import matplotlib.pyplot as plt
import math
from os import listdir
from os.path import isfile, join
import numpy as np
from numpy import pi
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

def mulliken_plot(filepath, filename="nah", energyshift=0, ymin=-5, ymax=5, substate=0, color_dict={}, eq=True, debug=True):
    ########################################
    # Data
    ########################################
    print("he")
    if debug:
        print("setting up workspace...")
    species_id = {}  # map species to its index in geometry.in: Pb --> 1
    sub_st = False
    if not substate == 0:
        sub_st = True  # 0, s; 1, p; 2, d; ...
    # line_energy_id = 1                                        # not currently used?
    line_mlk_start_id = 5  # this is the column where contributions are in mulliken outputs
    # line_atom_id = 3                                          # also not currently necessary

    ########################################
    # Figure settings
    ########################################
    matplotlib.rc('text', usetex='true')
    matplotlib.rcParams.update({'font.size': 20})
    matplotlib.rcParams['axes.linewidth'] = 2
    fig = plt.gcf()
    fig.set_size_inches(6, 6)
    linewidth = 1
    black_bands = 1
    pythoncolors = ['r', 'y', 'b', 'c', 'm', 'g', 'pink', 'purple', 'black']
    markersizeunit = 2
    ############################################
    # geometry.in
    ############################################
    if debug:
        print("reading data from ", filepath + "/geometry.in")
    atoms = []
    latvec = []
    rlatvec = []
    species_color = {}

    for line in open(filepath + "/geometry.in"):
        words = line.strip().split()
        if len(words) == 0:
            continue
        if words[0] == "lattice_vector":
            latvec.append([float(i) for i in words[1:4]])
        if line.strip().startswith("atom"):
            atoms.append(words[-1])  # full list of atoms in order they appear
    species = list(set(atoms))  # set of different types of atoms
    if color_dict == 0:
        for s in range(len(species)):
            species_color[species[s]] = pythoncolors[s]
    else:
        species_color = color_dict
    print("color key: ", species_color)

    for i in range(len(atoms)):  # both 2i and 2i+1 for the two spin states representing an atom
        if atoms[i] not in species_id:
            species_id[atoms[i]] = []
            species_id[atoms[i]].append(2 * i)
            species_id[atoms[i]].append(2 * i + 1)
        else:
            species_id[atoms[i]].append(2 * i)
            species_id[atoms[i]].append(2 * i + 1)

    # Calculate reciprocal lattice vectors
    volume = (np.dot(latvec[0], np.cross(latvec[1], latvec[2])))
    rlatvec.append(2 * pi * np.cross(latvec[1], latvec[2]) / volume)
    rlatvec.append(2 * pi * np.cross(latvec[2], latvec[0]) / volume)
    rlatvec.append(2 * pi * np.cross(latvec[0], latvec[1]) / volume)

    if debug:
        print("Lattice vectors:")
        for j in range(3):
            print(latvec[j])
        print("atoms:")
        print(atoms)
        print("species:")
        print(species)
        print("species_id:")
        print(species_id)
        print("species_color:", species_color)
        print("Reciprocal lattice vectors:")
        for j in range(3):
            print(rlatvec[j])
        print()

    ############################################
    # control.in
    ############################################
    if debug:
        print("reading data from ", filepath + "/control.in")
    else:
        print("Processing control.in and geometry.in...")
    kpoint = []
    band_len = []
    xvals = []
    band_len_tot = []

    for line in open(filepath + "/control.in"):
        #print(line.strip())
        if line.strip().startswith('output\tband'):
            words = line.strip().split()
            kpoint.append([float(i) for i in words[2:8]])
            n_sample = int(words[-3])  # n_sample is the number of integration points
    for i in kpoint:
        kvec = []
        xval = []
        for j in range(3):
            kvec.append(i[j + 3] - i[j])
        temp = math.sqrt(sum([k * k for k in list(np.dot(kvec, rlatvec))]))  # length of kpath segment
        step = temp / (n_sample - 1)
        for i in range(n_sample):
            xval.append(i * step)
        xvals.append(xval)  # list of all points calculated at
        band_len.append(temp)  # list of lengths of k-path segments
    all = sum(band_len)
    print("band lens: " + str(band_len))
    for i in range(len(band_len)):
        if i == 0:
            band_len_tot.append(0)
        else:
            band_len_tot.append(sum(band_len[:i]))
    for i in range(len(xvals)):
        xvals[i] = [j + band_len_tot[i] for j in
                    xvals[i]]  # basically separating each band path by assigning sequential locations on the x-axis
    if debug:
        print('K path')
        print(kpoint)
        print(n_sample)

    ############################################
    # band10**.out							   #
    ############################################
    print("Analyzing band100*.out files...")
    all_files = [join(filepath, f) for f in listdir(filepath) if isfile(join(filepath, f))]
    band_files = [f for f in all_files if 'band' in f and f.endswith('.out') and 'bandmlk' not in f]
    band_files.sort()
    band_mlkfiles = [f for f in all_files if 'bandmlk' in f and f.endswith('.out')]
    band_mlkfiles.sort()
    if len(band_mlkfiles):
        mlk = True
    else:
        mlk = False

    if debug:
        print(band_files)
        if mlk:
            print(band_mlkfiles)
            print("mlk: " + str(mlk))

    energys = []
    for file_id in range(len(band_files)):
        temp_en = []
        file = band_files[file_id]
        curTime = time.time()
        print("Processing " + file, end='\t\t')
        f = open(file)
        for line in f:
            words = line.strip().split()
            energy = []
            occ_ene = words[4:]
            # print occ_ene
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
        x_pts.append(all)
        # plt.yticks(range(ymin, ymax + 1), [])
        plt.xticks(ticks=x_pts, labels=('$X$', '$\Gamma$', '$Y\|L$', '$\Gamma$', '$K$', 'a', 'b'))
        plt.ylabel('Energy (eV)')
        plt.xlabel("Wave vector, $k$")
        # plt.title("3.9\% Doping", weight='bold')
        plt.title("(S-BrMBA)$_2$PbI$_4$", weight='bold')
        plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
        plt.tight_layout()
        plt.show()
        # plt.savefig(filename, dpi=300, bbox_inches='tight')
        return 0

    ############################################
    # mlk									   #
    ############################################
    print("Analyzing bandmlk100*.out files...")

    # for file_id in range(len(band_mlkfiles)):
    for file_id in range(len(band_mlkfiles)):
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
                        # if currentMlk[i] > 0:        # https://github.com/romankempt/aimstools/blob/c7ddf9b09c294f503687eac4e5a49e9cf9860469/aimstools/bandstructures/mulliken_bandstructure.py
                        # "Note: I'm forcing all l-contributions below zero to be zero, see discussion with Volker Blum."
                        spec_mlk[i_spec] += currentMlk[j]  # adding up all contributions of each atom type
                max = 0
                max_spec = 'Pb'
                for spec in species:
                    if spec_mlk[spec] > max:
                        max = spec_mlk[spec]
                        max_spec = spec
                if currentK % 4 == 0:
                    if max_spec == "Pb" or max_spec == "I":
                        plt.plot(xvals[file_id][currentK - 1], energys[file_id][currentK - 1][currentState - 1],
                                 species_color[max_spec] + 'o', markersize=markersizeunit,
                                 markeredgecolor=species_color[max_spec])
            i += 1
        print(str(time.time() - curTime) + " seconds")

    x_pts = []
    x_pts.append(0)
    for kpoint_x in band_len_tot[1:]:
        x_pts.append(kpoint_x)
        plt.axvline(kpoint_x, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
        plt.axhline(0, color='k', linestyle='--', lw=linewidth).set_dashes([5, 5])
    x_pts.append(all)
    # plt.yticks(range(ymin, ymax + 1), [])
    plt.xticks(ticks=x_pts, labels=('$X$', '$\Gamma$', '$Y\|L$', '$\Gamma$', '$K$'))
    plt.ylabel('Energy (eV)')
    plt.xlabel("Wave vector, $k$")
    # plt.title("3.9\% Doping", weight='bold')
    plt.title("(S-BrMBA)$_2$PbI$_4$", weight='bold')
    plt.axis([0, xvals[len(xvals) - 1][len(xvals[len(xvals) - 1]) - 1], ymin, ymax])
    plt.tight_layout()
    plt.show()
    # plt.savefig(filename, dpi=300, bbox_inches='tight')
