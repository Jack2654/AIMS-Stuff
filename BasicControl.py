# This file compiles basic functions for doing shit with control.in files
# In particular, the following functions exist here:
#
# -> write_control(geometry, base_control, defaults):
#                   given a geometry.in file, a template control.in, and a folder with species
#                   defaults, writes control.in file in same directory as geometry.in file
# -> write_all_controls(directory, base_control, defaults):
#                   using base_control as a template, writes a control file in each directory in the given directory
# -> species_default_path(element, defaults_folder):
#                   given an element and a path to the defaults, returns the path to the species default file
# -> list_of_defaults(filepath, defaul):
#                   returns a list of paths to all required defaults files for the given folder
# -> read_control_for_bands():
#                   reads control.in file as necessary for band outputs

# imports:
import os
import BasicGeo as bg
from PyAstronomy import pyasl
import math
import numpy as np


# functions:
def write_control(geometry, base_control, defaults, additional=""):
    if not os.path.exists(geometry):
        print("bad path given to write_control " + geometry)
        return 0
    defaults = list_of_defaults(geometry, defaults)
    with open(base_control, "r") as f:
        cont = f.read()
    with open(geometry[:-11] + "control.in", "w") as f:
        f.writelines(cont)
        f.write("\n" + additional + "\n")
        for default in defaults:
            with open(default, "r") as g:
                f.write(g.read())


def write_all_controls(directory, base_control, defaults, additional=""):
    for direct in os.listdir(os.fsencode(directory)):
        dir_name = os.fsdecode(direct)
        geometry = directory + dir_name + "/geometry.in"
        write_control(geometry, base_control, defaults, additional)


def species_default_path(element, defaults_folder):
    an = pyasl.AtomicNo()
    number = f"{an.getAtomicNo(element):02}"
    return defaults_folder + number + "_" + element + "_default"


def list_of_defaults(filepath, defaults_folder):
    defaults = []
    elements = bg.get_species_from_geo(filepath)
    for e in elements:
        defaults.append(species_default_path(e, defaults_folder))
    return defaults


def read_control_for_bands(filepath, rlatvec, bands, eq=False, debug=False):
    if debug:
        print("reading data from ", filepath + "/control.in")
    kpoint = []
    band_len = []
    xvals = []
    band_len_tot = []

    counter = 0
    if filepath[-1] == "/":
        read_path = filepath + "control.in"
    else:
        read_path = filepath + "/control.in"
    for line in open(read_path):
        if line.strip().startswith('output') and "band" in line:
            if counter in bands:
                words = line.strip().split()
                kpoint.append([float(i) for i in words[2:8]])
                n_sample = int(float((words[-3])))  # n_sample is the number of integration points
            counter += 1
    for i in kpoint:
        kvec = []
        xval = []
        for j in range(3):
            kvec.append(i[j + 3] - i[j])
        temp = math.sqrt(sum([k * k for k in list(np.dot(kvec, rlatvec))]))  # length of kpath segment
        if eq:
            step = 0.5 / (n_sample - 1)
        else:
            step = temp / (n_sample - 1)

        for i in range(n_sample):
            xval.append(i * step)
        xvals.append(xval)  # list of all points calculated at

        if eq:
            band_len.append(0.5)
        else:
            band_len.append(temp)  # list of lengths of k-path segments
    tot_len = sum(band_len)

    for i in range(len(band_len)):
        if i == 0:
            band_len_tot.append(0)
        else:
            band_len_tot.append(sum(band_len[:i]))
    for i in range(len(xvals)):
        xvals[i] = [j + band_len_tot[i] for j in
                    xvals[i]]  # basically separating each band path by assigning sequential locations on the x-axis
    if debug:
        print("band lens: " + str(band_len))
        print('K path')
        print(kpoint)
        print(n_sample)

    # xvals[0].reverse()
    # xvals[2].reverse()

    return xvals, band_len_tot, tot_len


def write_species(control_file, geo_file, species_folder):
    defaults = list_of_defaults(geo_file, species_folder)
    with open(control_file, "r") as f:
        ctrl = f.readlines()
    with open(control_file, "w") as f:
        f.writelines(ctrl)
        for spec in defaults:
            with open(spec, "r") as g:
                info = g.readlines()
            f.writelines(info)


def count_gaussian_bases(file):
    gauss_counts = {}
    gauss_count = -1
    nao_counts = {}
    nao_count = -1
    cur_spec = "X"
    include_nao = True

    with open(file, "r") as f:
        for line in f:
            if "include_min_basis" in line and ".false." in line:
                include_nao = False
            if "species" in line and "#" not in line:
                gauss_count = 0
                nao_count = 0
                cur_spec = line.split()[1]
                gauss_counts.update({cur_spec: 0})
                nao_counts.update({cur_spec: 0})
            if len(line.strip()) > 0:
                if "#" not in line.split()[0]:
                    if gauss_count > -1 and "gaussian" in line:
                        gauss_count += int(line.split()[2])
                        gauss_counts.update({cur_spec: gauss_count})
                    if nao_count > -1 and "valence" in line:
                        nao_count += 1
                        nao_counts.update({cur_spec: nao_count})
                    if nao_count > -1 and "hydro" in line:
                        nao_count += 1
                        nao_counts.update({cur_spec: nao_count})
                    if nao_count > -1 and "ionic" in line:
                        nao_count += 1
                        nao_counts.update({cur_spec: nao_count})

    res = []
    for key in gauss_counts.keys():
        if include_nao:
            temp = gauss_counts.get(key) + nao_counts.get(key)
        else:
            temp = gauss_counts.get(key)
        res.append(temp)
    return "\t".join([str(x) for x in res])


def recip_to_real(folder):
    control = folder + "control.in"
    geometry = folder + "geometry.in"
    rv = bg.reciprocal_vectors(geometry)
    bands = []
    with open(control, "r") as f:
        for line in f:
            if "band" in line:
                if len(line.split()) == 11:
                    temp = line.split()[2:8]
                    bands.append([float(x) for x in temp])
    for band in bands:
        direction = [0.0, 0.0, 0.0]
        for count, elem in enumerate(direction):
            direction[count] += rv[0][count] * band[0]
        print(direction)

