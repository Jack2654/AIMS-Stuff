# This file compiles basic functions for manipulating geometry.in files
# In particular, the following functions exist here:
#
# -> atoms(file):
#               returns an array with all given atoms in the file specified
# -> atoms_trimmed(file):
#               returns an array of arrays containing only the position information of each atom
# -> lattice_vectors(file):
#               returns an array of arrays of floats of the lattice vectors in the file specified
# -> distance(fileA, a, fileB, b):
#               returns the distance between points A and B in files A and B
# -> recenter(file, writepath, index, x, y, z):
#               re-centers a geometry.in file about the index provided using the coordinates provided
# -> reflect(file):
#               prints out some atoms reflected SHOULD BE REWORKED
# -> rot(p1, theta, axis):
#               literally no clue buster
# -> rotate(filepath, writepath, theta, axis):
#               rotates a geometry (about 0,0,0) about the given axis using the given theta
# -> move(read, write, dist):
#               moves a structure by an amount specified as (x,y,z) in dist
# -> xyz_to_geo(filepath, writepath):
#               translates a .xyz file to a .in file using ASE
# -> get_species_from_geo(filepath):
#               returns a set of the species present in a geometry.in file
#
#
#
#
#

# imports:
import math
from scipy.spatial.transform import Rotation as R
import VectorToolkit as vt
from ase.io import read, write


# functions:
def atoms(file):
    at = []
    with open(file, "r") as f:
        for line in f:
            if "atom" in line:
                at.append(line.strip())
    return at


def atoms_trimmed(file):
    at = atoms(file)
    ret = []
    for a in at:
        temp = a.split()
        ret.append([float(x) for x in temp[1:4]])
    return ret


def lattice_vectors(file):
    lv = []
    with open(file, "r") as f:
        for line in f:
            if "lattice" in line:
                lv.append(line)
    ret = []
    for lat in lv:
        temp = lat.split()
        ret.append([float(x) for x in temp[1:4]])
    return ret


def distance(fileA, a, fileB, b):
    x1 = atoms_trimmed(fileA)[a-1]
    x2 = atoms_trimmed(fileB)[b-1]
    return vt.dist(x1, x2)


def recenter(filepath, writepath, index, x, y, z):
    with open(filepath, "r") as f:
        j = 0
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
            else:
                st.append(ln)
        for i in lv:
            j += 1
            if j == index:
                temp = i.split()
                t1 = float(temp[0]) - x, float(temp[1]) - y, float(temp[2]) - z
                p1 = t1
    with open(writepath, "w") as f:
        f.writelines(st)
        for i in lv:
            temp = i.split()
            t1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = t1[0] - p1[0], t1[1] - p1[1], t1[2] - p1[2]
            f.write("atom " + str(final[0]) + " " + str(final[1]) + " " + str(final[2]) + " " + str(temp[3]) + "\n")


def reflect(filepath):
    at = atoms(filepath)
    for i in at:
        temp = i[5:].split()
        p1 = float(temp[0]), float(temp[1]), float(temp[2])
        final = p1[0], -p1[1], p1[2]
        print("atom", final[0], final[1], final[2], temp[3])


def rot(p1, theta, axis):
    rotation_radians = theta
    rotation_axis = axis
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    rotated_vec = rotation.apply(p1)
    return rotated_vec


def rotate(filepath, writepath, theta, axis):
    at = atoms(filepath)
    lv = lattice_vectors(filepath)
    with open(writepath, "w") as f:
        f.writelines(lv)
        for i in at:
            temp = i[5:].split()
            p1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = rot(p1, theta, axis)
            f.write("atom " + str(final[0]) + " " + str(final[1]) + " " + str(final[2]) + " " + str(temp[3]) + "\n")


def move(filepath, writepath, dist):
    lv = lattice_vectors(filepath)
    at_full = atoms(filepath)
    at = atoms_trimmed(filepath)
    ret = []
    for a in at:
        temp = a[0] + dist[0], a[1] + dist[1], a[2] + dist[2]
        ret.append(temp)
    with open(writepath, "w") as f:
        f.writelines(lv)
        i = 0
        for atm in ret:
            f.write("atom " + str(atm[0]) + " " + str(atm[1]) + " " + str(atm[2]) + " " + at_full[i].split()[4] + "\n")
            i += 1


def xyz_to_geo(filepath, writepath):
    a = read(filepath)
    a.write(writepath)


def get_species_from_geo(filepath):
    with open(filepath, "r") as f:
        species = []
        for ln in f:
            if "atom" in ln:
                species.append(ln.split()[4])
    return set(species)


def move_into_unit_cell(filepath, writepath):
    lines = []
    with open(filepath, "r") as f:
        for line in f:
            lines.append(line.strip())
    at_f = atoms(filepath)
    at = atoms_trimmed(filepath)
    lv = lattice_vectors(filepath)
    at_ret = []
    j = 0
    for a in at:
        i = 0
        at_ret.append([a[0], a[1], a[2]])
        for x in a:
            if x < 0:
                at_ret[j][0] += lv[i][0]
                at_ret[j][1] += lv[i][1]
                at_ret[j][2] += lv[i][2]
            if x > lv[i][i]:
                at_ret[j][0] -= lv[i][0]
                at_ret[j][1] -= lv[i][1]
                at_ret[j][2] -= lv[i][2]
            i += 1
        j += 1

    at_return = []
    i = 0
    for atom in at_ret:
        line = "atom " + str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + " " + at_f[i].split()[-1:][0] + "\n"
        at_return.append(line)
        i += 1

    with open(writepath, "w") as f:
        for ln in lines:
            if "atom" not in ln:
                f.write(ln + "\n")
        f.writelines(at_return)


