# This file compiles basic functions for manipulating geometry.in files
# In particular, the following functions exist here:
#
# -> atom(file, num):
#               returns an array specifying an atom from the given file
# -> atoms(file):
#               returns an array with all given atoms in the file specified
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
# ->
#
#
#
#
#
#
#
#
#
#

# imports:
import math
from scipy.spatial.transform import Rotation as R


def atom(file, num):
    at = []
    with open(file, "r") as f:
        for line in f:
            if "atom" in line:
                at.append(line)
    return at[num - 1]


def atoms(file):
    at = []
    with open(file, "r") as f:
        for line in f:
            if "atom" in line:
                at.append(line)
    return at


def distance(fileA, a, fileB, b):
    x1 = [float(x) for x in atom(fileA, a).split()[1:4]]
    x2 = [float(x) for x in atom(fileB, b).split()[1:4]]
    return math.sqrt(sum([(x1[x] - x2[x]) ** 2 for x in range(2)]))


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
                t1 = float(temp[0])-x, float(temp[1])-y, float(temp[2])-z
                p1 = t1
    with open(writepath, "w") as f:
        f.writelines(st)
        for i in lv:
            temp = i.split()
            t1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = t1[0]-p1[0],t1[1]-p1[1],t1[2]-p1[2]
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
    with open(filepath, "r") as f:
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
            elif ln.startswith("lattice"):
                st.append(ln)
    with open(writepath, "w") as f:
        f.writelines(st)
        for i in lv:
            temp = i.split()
            p1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = rot(p1, theta, axis)
            f.write("atom "+str(final[0])+" "+str(final[1])+" "+str(final[2])+" "+str(temp[3])+"\n")


