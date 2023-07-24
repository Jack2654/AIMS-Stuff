import matplotlib
import matplotlib.pyplot as plt
import math
from os import listdir
from os.path import isfile, join
import numpy as np
from numpy import pi
import time


# the band plotting functionality works, however mulliken output files generally do not contain all states calculated
# and thus the plot doesn't represent the full range of bands of a structure
def mull_plot(filepath=".", energy_shift=0, y_min=-5, y_max=5, color_dict={}, band_path=0, debug=False, quiet=False):
    line_width = 1
    black_bands = 1
    marker_size_unit = 2
    matplotlib.rc('text', usetex='true')
    matplotlib.rcParams.update({'font.size': 25})
    matplotlib.rcParams['axes.linewidth'] = 2
    fig = plt.gcf()
    fig.set_size_inches(6, 8)
    if debug:
        flags = 0       # full debug state, print all info
    elif not quiet:
        flags = 1       # normal run state, print out timing
    else:
        flags = 2       # quiet run state, print out nothing

    curTime = time.time()
    if flags == 0:
        print("Processing control.in and geometry.in...")
    elif flags == 1:
        print("Processing control.in and geometry.in...", end="\t\t")
    ############################################
    # geometry.in
    ############################################
    atoms = []
    lat_vec = []
    r_lat_vec = []
    species_id = {}
    for line in open(filepath + "/geometry.in"):
        words = line.strip().split()
        if len(words) == 0:
            continue
        if words[0] == "lattice_vector":
            lat_vec.append([float(i) for i in words[1:4]])
        if line.strip().startswith("atom"):
            atoms.append(words[-1])                     # full list of atoms in order they appear
    species = list(set(atoms))                          # set of different types of atoms
    for i in range(len(atoms)):                         # both 2i and 2i+1 for the two spin states representing an atom
        if atoms[i] not in species_id:
            species_id[atoms[i]] = []
        species_id[atoms[i]].append(2 * i)
        species_id[atoms[i]].append(2 * i + 1)
    volume = (np.dot(lat_vec[0], np.cross(lat_vec[1], lat_vec[2])))  # reciprocal lattice vector calculations
    r_lat_vec.append(2 * pi * np.cross(lat_vec[1], lat_vec[2]) / volume)
    r_lat_vec.append(2 * pi * np.cross(lat_vec[2], lat_vec[0]) / volume)
    r_lat_vec.append(2 * pi * np.cross(lat_vec[0], lat_vec[1]) / volume)
    ############################################
    # control.in
    ############################################
    k_point = []
    n_sample = []
    for line in open(filepath + "/control.in"):
        if line.strip().startswith('output band'):
            words = line.strip().split()
            k_point.append([float(i) for i in words[2:8]])
            n_sample.append(int(words[-3]))  # n_sample is the number of integration points
    if flags == 0:
        print("Lattice vectors:")
        for j in range(3):
            print(lat_vec[j])
        print("atoms:")
        print(atoms)
        print("species:")
        print(species)
        print("species_id:")
        print(species_id)
        print('K path')
        print(k_point)
        print('points along each path')
        print(n_sample)
        print("This took " + str(time.time() - curTime) + " seconds")
    elif flags == 1:
        print(str(time.time() - curTime) + " seconds")
    ############################################
    # variable setup
    ############################################
    species_color = {}
    band_len = []
    x_vals = []
    band_len_tot = []

    # this section finds bandmlk100*.out files and only selects the ones specified by band path to use
    all_files = [join(filepath, f) for f in listdir(filepath) if isfile(join(filepath, f))]
    band_mlk_files = [f for f in all_files if 'bandmlk' in f and f.endswith('.out')]
    band_mlk_files.sort()
    if band_path == 0:
        band_path = [i+1 for i in range(len(band_mlk_files))]
    band_mlk_files = [band_mlk_files[x] for x in range(len(band_mlk_files)) if (x+1) in band_path]
    k_point = [k_point[x] for x in range(len(k_point)) if (x+1) in band_path]

    for k, i in enumerate(k_point):
        k_vec = []
        x_val = []
        for j in range(3):
            k_vec.append(i[j + 3] - i[j])
        temp = math.sqrt(sum([k * k for k in list(np.dot(k_vec, r_lat_vec))]))  # length of kpath segment
        step = temp / (n_sample[k] - 1)
        for j in range(n_sample[k]):
            x_val.append(j * step)
        x_vals.append(x_val)  # list of all points calculated at
        band_len.append(temp)  # list of lengths of k-path segments
    for i in range(len(band_len)):
        if i == 0:
            band_len_tot.append(0)
        else:
            band_len_tot.append(sum(band_len[:i]))
    for i in range(len(x_vals)):
        x_vals[i] = [j + band_len_tot[i] for j in x_vals[i]]  # separating band paths by assigning sequential locations

    python_colors = ['r', 'y', 'b', 'c', 'm', 'g', 'pink', 'purple', 'black']
    if color_dict == 0:  # setting up colors for plotting
        for s in range(len(species)):
            species_color[species[s]] = python_colors[s]
    else:
        species_color = color_dict

    if flags == 0:
        print("color key: ", species_color)
        print("Reciprocal lattice vectors:")
        for j in range(3):
            print(r_lat_vec[j])
        print()
    ############################################
    # bandmlk100*.out
    ############################################
    line_mlk_start_id = 5  # within bandmlk100*.out files the 6th element in a line is the total angular momentum
    currentMlk = []
    for i, file in enumerate(band_mlk_files):
        curTime = time.time()
        all_states = []
        state_dict = {}
        state_added = False
        if flags == 0:
            print("Processing " + file)
        elif flags == 1:
            print("Processing " + file, end='\t\t')
        for line in open(file):
            if line.strip().startswith("k point number:"):
                currentK = int(line.strip().split()[3][:-1])
                kpt = line.strip()
            elif line.strip().startswith("State") and len(line.strip().split()) == 2:
                state_added = False
                words = line.strip().split()
                currentState = int(words[1])
                currentMlk = []
                if currentState not in state_dict:
                    state_dict[currentState] = []
                    all_states.append(currentState)
            elif line.strip() == '':
                continue
            elif line.strip()[0].isdigit():
                words = line.strip().split()
                if not state_added:
                    state_dict[currentState].append(float(words[1]))
                    state_added = True
                temp = 0                                   # note here, matching the usage of aimstools package, # all l < 0 are set to 0 then summed
                for x in words[line_mlk_start_id + 1:]:
                    if float(x) > 0:
                        temp += float(x)
                currentMlk.append(temp)
            if len(currentMlk) == len(atoms) * 2:
                continue
        for st in all_states:
            if 90 < i < 100: # not all states exist in mlk output files, idk which do :(
                print(i)
                plt.plot(x_vals[i], state_dict[st], color='k', lw=black_bands)
        print(str(time.time() - curTime) + " seconds")

    x_pts = [0]
    for k_point_x in band_len_tot[1:]:
        x_pts.append(k_point_x)
        plt.axvline(k_point_x, color='k', linestyle='--', lw=line_width).set_dashes([5, 5])
        plt.axhline(0, color='k', linestyle='--', lw=line_width).set_dashes([5, 5])
    x_pts.append(all)
    # plt.yticks(range(ymin, ymax + 1), [])
    plt.xticks(ticks=x_pts, labels=('$X$', '$\Gamma$', '$Y\|L$', '$\Gamma$', '$K$'))
    plt.ylabel('Energy (eV)')
    plt.xlabel("Wave vector, $k$")
    plt.title("3.9\% Doping")
    plt.axis([0, x_vals[len(x_vals) - 1][len(x_vals[len(x_vals) - 1]) - 1], y_min, y_max])
    plt.show()
    # plt.savefig(filename, dpi=300, bbox_inches='tight')