# This file compiles basic functions for manipulating geometry.in files
# In particular, the following functions exist here:
#
# -> atoms(file):
#               returns an array with all given atoms in the file specified
# -> atoms_trimmed(file):
#               returns an array of arrays containing only the position information of each atom
# -> lattice_vectors(file):
#               returns an array of arrays of floats of the lattice vectors in the file specified
# -> reciprocal_vectors(file):
#               returns an array of arrays of floats of the reciprocal lattice vectors based on the file specified
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
# -> move_into_unit_cell(filepath, writepath):
#               moves atoms outside the unit cell to inside the unit cell (approximately)
# -> read_geo_for_bands(filepath, color_dict={}, debug=False):
#               reads geometry.in file as needed for band outputs
# ->

# imports:
import math
import random
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from ase.io import read, write
import VectorToolkit as vt
import numpy as np
from numpy import pi
from mpl_toolkits.mplot3d import Axes3D


# functions:
def atoms(file):
    at = []
    with open(file, "r") as f:
        for line in f:
            if "atom" in line and "#" not in line:
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


def reciprocal_vectors(file):
    lv = lattice_vectors(file)
    a = lv[0]
    b = lv[1]
    c = lv[2]
    volume = abs(np.dot(a, np.cross(b, c)))
    B1 = tuple((2 * math.pi / volume) * x for x in np.cross(b, c))
    B2 = tuple((2 * math.pi / volume) * x for x in np.cross(c, a))
    B3 = tuple((2 * math.pi / volume) * x for x in np.cross(a, b))
    return [B1, B2, B3]


def distance(fileA, a, fileB, b):
    x1 = atoms_trimmed(fileA)[a - 1]
    x2 = atoms_trimmed(fileB)[b - 1]
    return vt.dist(x1, x2)


def recenter(filepath, writepath, index, x=0, y=0, z=0):
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
        lats = ["lattice_vector " + " ".join([str(y) for y in x]) + "\n" for x in lv]
        f.writelines(lats)
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
            if "atom" in ln and "#" not in ln:
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
    count = 0
    for a in at:
        i = 0
        at_ret.append([a[0], a[1], a[2]])
        for x in a:
            adj = False
            for lat in lv:
                if lat == lv[i]:
                    continue
                # if x > lat[i] + lv[i][i]:
                #    adj = True
                #    count += 1
                #    at_ret[j][0] -= lv[i][0]
                #    at_ret[j][1] -= lv[i][1]
                #    at_ret[j][2] -= lv[i][2]
            if x < 0 and not adj:
                count += 1
                at_ret[j][0] += lv[i][0]
                at_ret[j][1] += lv[i][1]
                at_ret[j][2] += lv[i][2]
            elif x > lv[i][i] and not adj:
                count += 1
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
    print(str(count) + " atoms were adjusted")


def read_geo_for_bands(filepath, color_dict={}, flags=1):
    if flags == 0:
        print("reading data from ", filepath + "/geometry.in")
    elif flags == 1:
        print("Processing control.in and geometry.in...")

    atoms = []
    latvec = []
    rlatvec = []
    species_id = {}  # map species to its index in geometry.in: Pb --> 1

    if filepath[-1] == "/":
        file = filepath + "geometry.in"
    else:
        file = filepath + "/geometry.in"
    if ".in" in filepath:
        file = filepath

    for line in open(file):
        words = line.strip().split()
        if len(words) == 0:
            continue
        if words[0] == "lattice_vector":
            latvec.append([float(i) for i in words[1:4]])
        if line.strip().startswith("atom"):
            atoms.append(words[-1])  # full list of atoms in order they appear
    species = list(set(atoms))  # set of different types of atoms

    for i in range(len(atoms)):  # both 2i and 2i+1 for the two spin states representing an atom
        if atoms[i] not in species_id:
            species_id[atoms[i]] = []
        species_id[atoms[i]].append(2 * i)
        species_id[atoms[i]].append(2 * i + 1)

    # Calculate reciprocal lattice vectors
    volume = (np.dot(latvec[0], np.cross(latvec[1], latvec[2])))
    rlatvec.append(2 * pi * np.cross(latvec[1], latvec[2]) / volume)
    rlatvec.append(2 * pi * np.cross(latvec[2], latvec[0]) / volume)
    rlatvec.append(2 * pi * np.cross(latvec[0], latvec[1]) / volume)

    if flags == 0:
        print("Lattice vectors:")
        for j in range(3):
            print(latvec[j])
        print("atoms:")
        print(atoms)
        print("species:")
        print(species)
        print("species_id:")
        print(species_id)
        print("species_color:", color_dict)
        print("Reciprocal lattice vectors:")
        for j in range(3):
            print(rlatvec[j])
        print()

    return rlatvec, atoms, species, species_id


def lattice_info(filename):
    lat = lattice_vectors(filename)
    mags = []
    for vec in lat:
        mag = math.sqrt((vec[0]) ** 2 + (vec[1]) ** 2 + (vec[2]) ** 2)
        mags.append(mag)
        # print(mag)
    # alpha:
    alpha = (180 / np.pi) * np.arccos(
        (lat[1][0] * lat[2][0] + lat[1][1] * lat[2][1] + lat[1][2] * lat[2][2]) / (mags[1] * mags[2]))
    print(alpha)
    # beta:
    beta = (180 / np.pi) * np.arccos(
        (lat[0][0] * lat[2][0] + lat[0][1] * lat[2][1] + lat[0][2] * lat[2][2]) / (mags[0] * mags[2]))
    print(beta)
    # gamma
    gamma = (180 / np.pi) * np.arccos(
        (lat[0][0] * lat[1][0] + lat[0][1] * lat[1][1] + lat[0][2] * lat[1][2]) / (mags[0] * mags[1]))
    print(gamma)


def delta_d_old(filename, center, edges, debug=False):
    lat = lattice_vectors(filename)
    at = atoms_trimmed(filename)
    at_center = at[center - 1]
    at_edges = [at[i - 1] for i in edges]
    at_edges_moved = []
    for atom in at_edges:
        for i in range(3):
            if atom[i] - at_center[i] > 5:
                atom = [atom[j] - lat[i][j] for j in range(3)]
            if atom[i] - at_center[i] < -5:
                atom = [atom[j] + lat[i][j] for j in range(3)]
        at_edges_moved.append(atom)

    distances = []
    for atom in at_edges_moved:
        distances.append(np.linalg.norm([atom[i] - at_center[i] for i in range(3)]))

    avg_distance = np.average(distances)
    result = 0
    for dist in distances:
        result += (dist - avg_distance) ** 2
    result *= (1 / 6) * (1 / (avg_distance ** 2))
    if debug:
        print('atom %.5f %.5f %.5f M' % (at_center[0], at_center[1], at_center[2]))
        for atom in at_edges_moved:
            print('atom %.5f %.5f %.5f H' % (atom[0], atom[1], atom[2]))
    return result


def sigma_squared_old(filename, center, edges, debug=False):
    # order of edges is [left, up, right, down, above, below]
    lat = lattice_vectors(filename)
    at = atoms_trimmed(filename)
    at_center = at[center - 1]
    at_edges = [at[i - 1] for i in edges]
    at_edges_moved = []
    for atom in at_edges:
        for i in range(3):
            if atom[i] - at_center[i] > 5:
                atom = [atom[j] - lat[i][j] for j in range(3)]
            if atom[i] - at_center[i] < -5:
                atom = [atom[j] + lat[i][j] for j in range(3)]
        at_edges_moved.append(atom)

    angles = [angle(at_edges_moved[0], at_center, at_edges_moved[1]),
              angle(at_edges_moved[0], at_center, at_edges_moved[4]),
              angle(at_edges_moved[0], at_center, at_edges_moved[5]),
              angle(at_edges_moved[0], at_center, at_edges_moved[3]),
              angle(at_edges_moved[1], at_center, at_edges_moved[4]),
              angle(at_edges_moved[1], at_center, at_edges_moved[5]),
              angle(at_edges_moved[1], at_center, at_edges_moved[2]),
              angle(at_edges_moved[2], at_center, at_edges_moved[4]),
              angle(at_edges_moved[2], at_center, at_edges_moved[5]),
              angle(at_edges_moved[2], at_center, at_edges_moved[3]),
              angle(at_edges_moved[3], at_center, at_edges_moved[4]),
              angle(at_edges_moved[3], at_center, at_edges_moved[5])]
    # print(angles)
    result = 0
    for ang in angles:
        result += (ang - 90) ** 2
    result /= 12

    return result


# returns angle in degrees
# can work for nD points
def angle(pt1, center, pt2):
    a = np.array(pt1)
    b = np.array(center)
    c = np.array(pt2)
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    result = np.arccos(cosine_angle)
    return np.degrees(result)


def angle_info(readfile, pts, shiftmap, direction=0):
    debug = False
    # pulls information from geometry.in and pulls three desired points
    lv = lattice_vectors(readfile)
    at = atoms_trimmed(readfile)
    points = [at[pts[0] - 1], at[pts[1] - 1], at[pts[2] - 1]]

    # moves points so they are positioned adequately
    if not shiftmap == -1:
        for i in range(3):
            for j in range(3):
                points[i][0] += shiftmap[i][j] * lv[j][0]
                points[i][1] += shiftmap[i][j] * lv[j][1]
                points[i][2] += shiftmap[i][j] * lv[j][2]

    p1, p2, p3 = points
    # computes normal vector to the "in-plane" plane
    connect_vect = [p3[i] - p1[i] for i in range(3)]
    b = math.sqrt((connect_vect[0] ** 2) / (connect_vect[0] ** 2 + connect_vect[1] ** 2))
    a = math.sqrt(1 - b ** 2)
    horizontal_vec = [a, b, 0]
    normal_vec = np.cross(connect_vect, horizontal_vec)
    normal_vec = [normal_vec[i] / np.linalg.norm(normal_vec) for i in range(3)]

    # computes projection of p2 onto "in-plane" plane
    p2_from_origin = [p2[i] - p1[i] for i in range(3)]
    in_plane_proj_factor = np.dot(normal_vec, p2_from_origin)
    in_plane_proj_vector = [in_plane_proj_factor * normal_vec[i] for i in range(3)]
    in_plane_point = [p2[i] - in_plane_proj_vector[i] for i in range(3)]

    # computes projection of p2 onto "out-of-plane" plane
    horizontal_vec = [horizontal_vec[i] / np.linalg.norm(horizontal_vec) for i in range(3)]
    p2_from_origin = [p2[i] - p1[i] for i in range(3)]
    out_of_plane_proj_factor = np.dot(horizontal_vec, p2_from_origin)
    out_of_plane_proj_vector = [out_of_plane_proj_factor * horizontal_vec[i] for i in range(3)]
    out_of_plane_point = [p2[i] - out_of_plane_proj_vector[i] for i in range(3)]

    if debug:
        print("Point 1: " + str(p1))
        print("Point 2: " + str(p2))
        print("Point 3: " + str(p3))
        print("In plane projection: " + str(in_plane_point))
        print("Out of plane projection: " + str(out_of_plane_point))
        print('Full Angle: %.10f' % angle(p1, p2, p3))
        print('In-plane angle: %.10f' % angle(p1, in_plane_point, p3))
        print('Out-of-plane angle: %.10f' % angle(p1, out_of_plane_point, p3))
    else:
        res = []
        if direction == 0:
            res.append(angle(p1, p2, p3))
            return res

        d = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
        if direction < 0:
            d *= -1
        if d > 0:
            res.append(angle(p1, p2, p3))
        elif d < 0:
            res.append(360 - angle(p1, p2, p3))
        else:
            res.append(0)

        d = (in_plane_point[0] - p1[0]) * (p3[1] - p1[1]) - (in_plane_point[1] - p1[1]) * (p3[0] - p1[0])
        if direction < 0:
            d *= -1
        if d > 0:
            res.append(angle(p1, in_plane_point, p3))
        elif d < 0:
            res.append(360 - angle(p1, in_plane_point, p3))
        else:
            res.append(0)

        d = (out_of_plane_point[0] - p1[0]) * (p3[1] - p1[1]) - (out_of_plane_point[1] - p1[1]) * (p3[0] - p1[0])
        if direction < 0:
            d *= -1
        if d > 0:
            res.append(angle(p1, out_of_plane_point, p3))
        elif d < 0:
            res.append(360 - angle(p1, out_of_plane_point, p3))
        else:
            res.append(0)

        return res


def angle_info_old(readfile, pts):
    debug = False
    i1 = pts[0]
    i2 = pts[1]
    i3 = pts[2]
    # pulls information from geometry.in and pulls three desired points
    lv = lattice_vectors(readfile)
    at = atoms_trimmed(readfile)
    p1 = at[i1 - 1]
    p2 = at[i2 - 1]
    p3 = at[i3 - 1]

    # moves points so they are positioned adequately

    # computes normal vector to the "in-plane" plane
    connect_vect = [p3[i] - p1[i] for i in range(3)]
    b = math.sqrt((connect_vect[0] ** 2) / (connect_vect[0] ** 2 + connect_vect[1] ** 2))
    a = math.sqrt(1 - b ** 2)
    horizontal_vec = [a, b, 0]
    normal_vec = np.cross(connect_vect, horizontal_vec)
    normal_vec = [normal_vec[i] / np.linalg.norm(normal_vec) for i in range(3)]

    # computes projection of p2 onto "in-plane" plane
    p2_from_origin = [p2[i] - p1[i] for i in range(3)]
    in_plane_proj_factor = np.dot(normal_vec, p2_from_origin)
    in_plane_proj_vector = [in_plane_proj_factor * normal_vec[i] for i in range(3)]
    in_plane_point = [p2[i] - in_plane_proj_vector[i] for i in range(3)]

    # computes projection of p2 onto "out-of-plane" plane
    horizontal_vec = [horizontal_vec[i] / np.linalg.norm(horizontal_vec) for i in range(3)]
    p2_from_origin = [p2[i] - p1[i] for i in range(3)]
    out_of_plane_proj_factor = np.dot(horizontal_vec, p2_from_origin)
    out_of_plane_proj_vector = [out_of_plane_proj_factor * horizontal_vec[i] for i in range(3)]
    out_of_plane_point = [p2[i] - out_of_plane_proj_vector[i] for i in range(3)]

    if debug:
        print("Point 1: " + str(p1))
        print("Point 2: " + str(p2))
        print("Point 3: " + str(p3))
        print("In plane projection: " + str(in_plane_point))
        print("Out of plane projection: " + str(out_of_plane_point))
        print('Full Angle: %.10f' % angle(p1, p2, p3))
        print('In-plane angle: %.10f' % angle(p1, in_plane_point, p3))
        print('Out-of-plane angle: %.10f' % angle(p1, out_of_plane_point, p3))
    else:
        res = []
        res.append(angle(p1, p2, p3))
        res.append(angle(p1, in_plane_point, p3))
        res.append(angle(p1, out_of_plane_point, p3))
        return res


def disturb_positions(readfile, writefile, max_disturbance, change_Cs=False):
    lv = lattice_vectors(readfile)
    at = atoms(readfile)
    with open(writefile, "w") as f:
        for lat in lv:
            f.write('lattice_vector\t%.8f\t%.8f\t%.8f\n' % (lat[0], lat[1], lat[2]))
        for atom in at:
            temp = atom.split()
            if not change_Cs and "Cs" in temp:
                # if -0.1 < float(temp[3]) < 0.1: <- appears to be wrong
                f.write('atom\t%s\t%s\t%s\t%s\n' % (temp[1], temp[2], temp[3], temp[-1]))
                continue
            coords = [float(x) for x in temp[1:4]]
            for i in range(3):
                coords[i] += max_disturbance * ((2 * random.random()) - 1)
            f.write('atom\t%.8f\t%.8f\t%.8f\t%s\n' % (coords[0], coords[1], coords[2], temp[-1]))


def maurer_displacement(readfile, corners, shiftmap, center, Ag=True):
    atoms = atoms_trimmed(readfile)
    lv = lattice_vectors(readfile)
    corner_loc = []
    for i in range(len(corners)):
        temp = atoms[corners[i] - 1]
        for j in range(3):
            temp[0] += shiftmap[i][j] * lv[j][0]
            temp[1] += shiftmap[i][j] * lv[j][1]
            temp[2] += shiftmap[i][j] * lv[j][2]
        corner_loc.append(temp)
    centroid = [0, 0, 0]
    for corner in corner_loc:
        centroid[0] += corner[0]
        centroid[1] += corner[1]
        centroid[2] += corner[2]
    centroid[0] /= len(corners)
    centroid[1] /= len(corners)
    centroid[2] /= len(corners)
    center_actual = atoms[center - 1]
    displacement = [centroid[i] - center_actual[i] for i in range(3)]
    displacement[2] = 0

    if Ag:
        index_a = 0
        index_b = 1
    else:
        index_a = 2
        index_b = 3
    arm_one = [corner_loc[index_a][x] - center_actual[x] for x in range(3)]
    arm_two = [corner_loc[index_b][x] - center_actual[x] for x in range(3)]
    diag_one = [arm_one[x] + arm_two[x] for x in range(3)]
    diag_two = [arm_two[x] - arm_one[x] for x in range(3)]
    diag_one = [x / np.linalg.norm(diag_one) for x in diag_one]
    diag_two = [x / np.linalg.norm(diag_two) for x in diag_two]

    # print(f'Centroid:\t\t%0.5f %0.5f %0.5f' % (centroid[0], centroid[1], centroid[2]))
    # print(f'Center:\t\t\t%0.5f %0.5f %0.5f' % (center_actual[0], center_actual[1], center_actual[2]))
    # print(f'Displacement:\t%0.5f %0.5f %0.5f' % (displacement[0], displacement[1], displacement[2]))
    # print(f'Diag One:\t%0.5f %0.5f %0.5f' % (diag_one[0], diag_one[1], diag_one[2]))
    # print(f'Diag Two:\t%0.5f %0.5f %0.5f\n' % (diag_two[0], diag_two[1], diag_two[2]))

    proj_one = abs(np.dot(displacement, diag_one))
    proj_two = abs(np.dot(displacement, diag_two))

    # code to verify projections produce original result, remove abs()
    # diag_one = [x * proj_one for x in diag_one]
    # diag_two = [x * proj_two for x in diag_two]
    # total = [(diag_one[i] + diag_two[i]) * np.linalg.norm(displacement) for i in range(3)]
    # print(total)
    # return f'%0.8f %0.8f' % (proj_one, proj_two)

    return proj_one, proj_two


def sigma_squared(readfile, corners, shiftmap, center):
    # order of corners is left, down, right, up, below, above
    atoms = atoms_trimmed(readfile)
    lv = lattice_vectors(readfile)
    corner_loc = []
    for i in range(len(corners)):
        temp = atoms[corners[i] - 1]
        for j in range(3):
            temp[0] += shiftmap[i][j] * lv[j][0]
            temp[1] += shiftmap[i][j] * lv[j][1]
            temp[2] += shiftmap[i][j] * lv[j][2]
        corner_loc.append(temp)
    center_actual = atoms[center - 1]

    angles = [angle(corner_loc[0], center_actual, corner_loc[1]),
              angle(corner_loc[0], center_actual, corner_loc[3]),
              angle(corner_loc[0], center_actual, corner_loc[4]),
              angle(corner_loc[0], center_actual, corner_loc[5]),
              angle(corner_loc[1], center_actual, corner_loc[2]),
              angle(corner_loc[1], center_actual, corner_loc[4]),
              angle(corner_loc[1], center_actual, corner_loc[5]),
              angle(corner_loc[2], center_actual, corner_loc[3]),
              angle(corner_loc[2], center_actual, corner_loc[4]),
              angle(corner_loc[2], center_actual, corner_loc[5]),
              angle(corner_loc[3], center_actual, corner_loc[4]),
              angle(corner_loc[3], center_actual, corner_loc[5])]
    result = 0
    for ang in angles:
        result += (ang - 90) ** 2
    result /= 12

    return result


def delta_d(readfile, corners, shiftmap, center):
    # order of corners is left, down, right, up, below, above
    atoms = atoms_trimmed(readfile)
    lv = lattice_vectors(readfile)
    corner_loc = []
    for i in range(len(corners)):
        temp = atoms[corners[i] - 1]
        for j in range(3):
            temp[0] += shiftmap[i][j] * lv[j][0]
            temp[1] += shiftmap[i][j] * lv[j][1]
            temp[2] += shiftmap[i][j] * lv[j][2]
        corner_loc.append(temp)
    center_actual = atoms[center - 1]

    distances = []
    for atom in corner_loc:
        distances.append(np.linalg.norm([atom[i] - center_actual[i] for i in range(3)]))

    avg_distance = np.average(distances)
    result = 0
    for dist in distances:
        result += (dist - avg_distance) ** 2
    result *= (1 / 6) * (1 / (avg_distance ** 2))
    return result


def geo_to_poscar(readfile, writefile):
    at = atoms(readfile)
    lv = lattice_vectors(readfile)
    with open(writefile, "w") as f:
        f.write(readfile + "\n")  # header comment
        f.write("1.00\n")  # scaling factor
        for lat in lv:
            f.write("\t".join(str(x) for x in lat) + "\n")

        atom_counts = {}
        for a in at:
            temp = a.split()[-1]
            if temp in atom_counts.keys():
                atom_counts[temp] = atom_counts[temp] + 1
            else:
                atom_counts[temp] = 1

        f.write("\t".join(atom_counts.keys()) + "\n")
        f.write("\t".join(str(x) for x in [atom_counts[y] for y in atom_counts.keys()]) + "\n")
        f.write("Cartesian\n")

        for key in atom_counts.keys():
            for a in at:
                if a.split()[-1] == key:
                    f.write("\t".join(str(x) for x in a.split()[1:4]) + "\n")

    with open(writefile, "r") as f:
        lines = f.readlines()
    for ln in lines:
        print(ln.strip())


def closest(corner, corners):
    min_index = -1
    min_value = 200
    second_min_index = -1
    second_min_value = 200
    for count, c in enumerate(corners):
        temp_dist = vt.dist(corner, c)
        if temp_dist != 0:
            if temp_dist < second_min_value:
                second_min_value = temp_dist
                second_min_index = count
            if temp_dist < min_value:
                second_min_value = min_value
                second_min_index = min_index
                min_value = temp_dist
                min_index = count
    return [min_index, second_min_index]


# writes new geometry file to desired location where the list of atoms have their species changed to the changed field
def rename_species(readfile, writefile, indices, original, changed):
    lv = lattice_vectors(readfile)
    full_lv = ["lattice_vector " + "\t".join([str(y) for y in x]) + "\n" for x in lv]
    at = atoms(readfile)
    at_expanded = [x.split() for x in at]
    for i, atom in enumerate(at_expanded):
        if i + 1 in indices:
            if atom[-1] == original:
                atom[-1] = changed
            else:
                print("CHECK YOUR INDEXES!!!")
    with open(writefile, "w") as f:
        f.writelines(full_lv)
        f.writelines(["\t".join(x) + "\n" for x in at_expanded])


def random_angle_info(readfile, method="normal", bonds=""):
    if bonds == "":
        point_sets = [(4, 7, 1), (4, 10, 1), (4, 8, 1), (4, 9, 1)]
        shiftmaps = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [0, 1, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 1, 0]]]
    else:
        point_sets = [[atom[0] for atom in bond] for bond in bonds]
        shiftmaps = [[[0, 0, 0] if len(atom) == 1 else [int(x) for x in list(atom[1])] for atom in bond] for bond in
                     bonds]
    angles = []
    # direction 1 = vertical, 2 = horizontal
    directions = [1, 2, 2, 1]
    # 1 = should be + coordinate direction, -1 = should be - coordinate direction
    desired_polarity = [1, 1, -1, -1]
    for count, points in enumerate(point_sets):
        debug = False
        # pulls information from geometry.in and pulls three desired points
        lv = lattice_vectors(readfile)
        at = atoms_trimmed(readfile)
        pts = [at[points[0] - 1], at[points[1] - 1], at[points[2] - 1]]
        # moves points so they are positioned adequately
        for i in range(3):
            for j in range(3):
                pts[i][0] += shiftmaps[count][i][j] * lv[j][0]
                pts[i][1] += shiftmaps[count][i][j] * lv[j][1]
                pts[i][2] += shiftmaps[count][i][j] * lv[j][2]

        if method == "flat":
            pts[0][2] = 0
            pts[1][2] = 0
            pts[2][2] = 0

        p1, p2, p3 = pts

        if method == "ignore":
            angles.append(angle(p1, p2, p3))
            continue

        slope = (p3[1] - p1[1]) / (p3[0] - p1[0])
        intercept = p3[1] - p3[0] * slope

        # case of vertical angle
        if directions[count] == 1:
            polarity = 1 if p2[0] > (p2[1] - intercept) / slope else -1
            temp_angle = angle(p1, p2, p3) if polarity == desired_polarity[count] else 360 - angle(p1, p2, p3)
        # case of vertical angle
        elif directions[count] == 2:
            polarity = 1 if p2[1] > p2[0] * slope + intercept else -1
            temp_angle = angle(p1, p2, p3) if polarity == desired_polarity[count] else 360 - angle(p1, p2, p3)
        angles.append(temp_angle)
        # print(temp_angle)
        # x_vals = [p1[0], p3[0]]
        # plt.plot(x_vals, [x * slope + intercept for x in x_vals], label=temp_angle)
        # plt.scatter([p1[0], p3[0]], [p1[1], p3[1]], color='green', label="_")
        # if polarity == desired_polarity[count]:
        #     plt.scatter(p2[0], p2[1], color='b', label="_")
        # else:
        #     plt.scatter(p2[0], p2[1], color='r', label="_")
    # plt.xlim(-7, 7)
    # plt.ylim(-1, 13)
    # plt.legend()
    # plt.show()
    return angles


def flatten(readfile, writefile):
    lv = lattice_vectors(readfile)
    at = atoms(readfile)
    with open(writefile, "w") as f:
        for lat in lv:
            f.write("lattice_vector " + " ".join([str(x) for x in lat]) + "\n")
        for atom in at:
            temp = atom.split()
            value = abs(float(temp[3]))
            if value < 2:
                f.write(" ".join(temp[:3]) + " 0.0 " + temp[-1] + "\n")
            else:
                f.write(atom + "\n")


# computes the four components needed for the Delta_L descriptor as described in the Group Meeting v4 Slides
def Delta_L(readfile, mode="Delta", bonds="", flag=""):
    if bonds == "":
        point_sets = [(4, 7, 1), (4, 10, 1), (4, 8, 1), (4, 9, 1)]
        shiftmaps = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [0, 1, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 1, 0]]]
    else:
        point_sets = [[atom[0] for atom in bond] for bond in bonds]
        shiftmaps = [[[0, 0, 0] if len(atom) == 1 else [int(x) for x in list(atom[1])] for atom in bond] for bond in
                     bonds]

    L_vals = []
    lengths = []
    lv = lattice_vectors(readfile)
    for count, points in enumerate(point_sets):
        at = atoms_trimmed(readfile)
        pts = [at[points[0] - 1], at[points[1] - 1], at[points[2] - 1]]
        # moves points so they are positioned adequately
        for i in range(3):
            for j in range(3):
                pts[i][0] += shiftmaps[count][i][j] * lv[j][0]
                pts[i][1] += shiftmaps[count][i][j] * lv[j][1]
                pts[i][2] += shiftmaps[count][i][j] * lv[j][2]
        p1, p2, p3 = pts
        L_1 = np.linalg.norm([p1[i] - p2[i] for i in range(3)])
        L_2 = np.linalg.norm([p2[i] - p3[i] for i in range(3)])
        if flag == "debug":
            print(pts)
        if mode == "Delta":
            L_vals.append((L_1 - L_2) ** 2)
        elif mode == "L2":
            L_vals.append((L_1 - 2.88196731) ** 2)
            L_vals.append((L_2 - 2.88196731) ** 2)
            lengths.append(L_1)
            lengths.append(L_2)
    if mode == "L2" and not bonds == "":
        avg_len = np.average(lengths)
        L_vals = []
        for leng in lengths:
            L_vals.append((leng - avg_len) ** 2)
    return L_vals


# computes the L_diff descriptor for AgBi structures
def L_diff(readfile, bonds="", mode="", frac=False):
    if bonds == "":
        point_sets = [(4, 7, 1), (4, 10, 1), (4, 8, 1), (4, 9, 1)]
        shiftmaps = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 0, 0]],
                     [[0, 0, 0], [0, 0, 0], [0, 1, 0]],
                     [[0, 0, 0], [0, 0, 0], [1, 1, 0]]]
    else:
        point_sets = [[atom[0] for atom in bond] for bond in bonds]
        shiftmaps = [[[0, 0, 0] if len(atom) == 1 else [int(x) for x in list(atom[1])] for atom in bond] for bond in
                     bonds]

    Ag_lens = []
    Bi_lens = []
    lv = lattice_vectors(readfile)
    for count, points in enumerate(point_sets):
        at = atoms_trimmed(readfile)
        pts = [at[points[0] - 1], at[points[1] - 1], at[points[2] - 1]]
        if frac:
            for i in range(3):
                for j in range(3):
                    pts[i][j] = pts[i][0] * lv[0][j] + pts[i][1] * lv[1][j] + pts[i][2] * lv[2][j]

        # moves points so they are positioned adequately
        for i in range(3):
            for j in range(3):
                pts[i][0] += shiftmaps[count][i][j] * lv[j][0]
                pts[i][1] += shiftmaps[count][i][j] * lv[j][1]
                pts[i][2] += shiftmaps[count][i][j] * lv[j][2]
        p1, p2, p3 = pts
        Ag_lens.append(np.linalg.norm([p1[i] - p2[i] for i in range(3)]))
        Bi_lens.append(np.linalg.norm([p2[i] - p3[i] for i in range(3)]))

    if mode == "all":
        avg_len = np.average(Ag_lens + Bi_lens)
        result = np.average([abs((elem - avg_len) / avg_len) for elem in Ag_lens + Bi_lens])
        return result

    Ag_avg = np.average(Ag_lens)
    Bi_avg = np.average(Bi_lens)
    Ag_diff = sum([abs((Ag_l - Ag_avg) / Ag_avg) for Ag_l in Ag_lens])
    Bi_diff = sum([abs((Bi_l - Bi_avg) / Bi_avg) for Bi_l in Bi_lens])
    result = (Ag_diff + Bi_diff) / 8
    return result


# computes the average in-plane Ag-halide length for AgBi-I model structures
def Ag_len(readfile):
    # want 4 -> 7, 8, 9, 10
    # in base structure length is 2.890625
    Ag_lens = []
    lv = lattice_vectors(readfile)
    at = atoms_trimmed(readfile)
    center = at[4 - 1]
    for point in [7, 8, 9, 10]:
        current = at[point - 1]
        Ag_lens.append((np.linalg.norm([center[i] - current[i] for i in range(3)]) - 2.890625) ** 2)
    return sum(Ag_lens)


def dict_to_in(input_dict):
    for atom in input_dict.keys():
        temp = input_dict[atom]
        print(f'atom {temp[0]} {temp[1]} {temp[2]} Bi')


def pull_data(readfile, bonds, shiftmap):
    # pulls lattice vectors and atoms from the geometry.in file
    lv = lattice_vectors(readfile)
    at = atoms_trimmed(readfile)
    atom_dict = {}
    for bond in bonds:
        for atom in bond:
            atom_dict[atom] = at[int(atom.split("_")[0]) - 1][:] if "_" in atom else at[int(atom) - 1][:]

    # shifts the atoms as necessary
    if shiftmap is not None:
        for shift_index in shiftmap.keys():
            current = atom_dict[shift_index]
            shifts = [int(x) for x in shiftmap[shift_index].split(",")]
            # iterate through x,y,z coordinates (i) and for each iterate through lattice vectors (j)
            for i in range(3):
                for j in range(3):
                    current[i] += shifts[j] * lv[j][i]
            atom_dict[shift_index] = current

    # recenters the atoms about the first B-site
    recenter = np.array(atom_dict[bonds[0][0]])
    for atom in atom_dict.keys():
        atom_dict[atom] = np.array(atom_dict[atom]) - recenter

    return atom_dict

    # rotates the atoms so the third B-site is also at z=0
    axis = atom_dict[bonds[0][2]] - atom_dict[bonds[2][2]]
    axis = axis / np.linalg.norm(axis)
    cur_flat_dist = math.sqrt(atom_dict[bonds[1][2]][0] ** 2 + atom_dict[bonds[1][2]][1] ** 2)
    cur_vert_dist = atom_dict[bonds[1][2]][2]
    cur_angle = np.arctan(cur_vert_dist / cur_flat_dist)
    rotation = R.from_rotvec(-1 * cur_angle * axis)  # -1 to "undo" the current rotation

    # Apply the rotation to all points
    for atom in atom_dict.keys():
        atom_dict[atom] = rotation.apply(atom_dict[atom])

    # rotates the atoms to minimize the height of the remaining 2 B-sites
    axis = atom_dict[bonds[1][2]] - atom_dict[bonds[0][0]]
    axis = axis / np.linalg.norm(axis)
    cur_flat_dist_1 = math.sqrt(atom_dict[bonds[1][2]][0] ** 2 + atom_dict[bonds[1][2]][1] ** 2)
    cur_vert_dist_1 = atom_dict[bonds[1][2]][2]
    cur_flat_dist_2 = math.sqrt(atom_dict[bonds[1][2]][0] ** 2 + atom_dict[bonds[1][2]][1] ** 2)
    cur_vert_dist_2 = atom_dict[bonds[1][2]][2]
    cur_angle = np.arctan(cur_vert_dist / cur_flat_dist)

    # not completed
    return atom_dict


def find_db_in(bond_coords):
    # first flatten all angles
    flat_coords = [[np.array([atom[0], atom[1], 0]) for atom in bond] for bond in bond_coords]
    # next compute the 'centroid' of metal B-sites
    centroid = np.array([0.0, 0.0, 0.0])
    for i in range(4):
        centroid += bond_coords[i][0]
    centroid = [centroid[i] / 4 for i in range(len(centroid))]
    plt.scatter(centroid[0], centroid[1], color='silver')
    beta_in = []
    # deal with each angle separately
    for count, bond in enumerate(flat_coords):
        # compute equation connecting metal B-sites
        m = (bond[0][1] - bond[2][1]) / (bond[0][0] - bond[2][0])
        b = bond[0][1] - m * bond[0][0]
        plt.scatter([bond[0][0], bond[2][0]], [bond[0][1], bond[2][1]], color='k')
        plt.plot([bond[0][0], bond[2][0]], [m * bond[0][0] + b, m * bond[2][0] + b], color='k')
        # project the halide onto the segment
        x_proj = (bond[1][0] + m * bond[1][1] - m * b) / (1 + m * m)
        y_proj = m * x_proj + b
        dist_actual = np.linalg.norm(bond[1] - centroid)
        dist_proj = np.linalg.norm(np.array([x_proj, y_proj, 0]) - centroid)
        # frame of reference has even bonds (0 and 2) inside and
        # odd bonds (1 and 3) outside
        cur_angle = angle(*bond)
        if dist_actual <= dist_proj:
            plt.scatter(bond[1][0], bond[1][1], color='green')
            beta_in.append(cur_angle if (count % 2) == 0 else 360 - cur_angle)
        else:
            plt.scatter(bond[1][0], bond[1][1], color='red')
            beta_in.append(cur_angle if (count % 2) == 1 else 360 - cur_angle)
    # plt.show()
    return beta_in


def find_db_out(bond_coords):
    # construct a plane connecting endpoints and parallel to (0, 0, 1)
    def construct_plane_normal(p1, p2):
        v1 = np.array(p2) - np.array(p1)
        v2 = np.array([0, 0, 1])
        n = np.cross(v1, v2)
        return n

    # project a point onto a given plane
    def project_point_to_plane(p, p0, n):
        w = np.array(p) - np.array(p0)
        scalar = np.dot(n, w) / np.dot(n, n)
        proj = np.array(p) - scalar * n
        return proj

    # project the point on the plane onto the line connecting endpoints, this helps with the frame of reference
    def project_point_onto_line(p1, p2, center):
        d = p1 - p2  # Direction vector
        AP = center - p1  # Vector from A to P
        time = np.dot(AP, d) / np.dot(d, d)  # Scalar projection
        proj = p1 + time * d  # Projected point
        return proj

    beta_out = []
    polarity = [0, 0, 0, 0]  # 1 = above, -1 = below
    for count, bond in enumerate(bond_coords):
        p1 = bond[0]
        p2 = bond[2]
        p_to_project = bond[1]
        normal = construct_plane_normal(p1, p2)
        p_proj = project_point_to_plane(p_to_project, p1, normal)
        beta_out.append(angle(p1, p_proj, p2))
        polarity[count] = 1 if p_proj[2] > project_point_onto_line(p1, p2, p_proj)[2] else -1
    untouched_beta = beta_out[:]

    if polarity[0] == polarity[1]:
        beta_out[2] = beta_out[2] if polarity[2] != polarity[0] else 360 - beta_out[2]
        beta_out[3] = beta_out[3] if polarity[3] != polarity[0] else 360 - beta_out[3]
        db_out = (beta_out[0] + beta_out[3] - beta_out[2] - beta_out[1]) / 2
    elif polarity[1] == polarity[2]:
        beta_out[0] = beta_out[0] if polarity[0] != polarity[1] else 360 - beta_out[0]
        beta_out[3] = beta_out[3] if polarity[3] != polarity[1] else 360 - beta_out[3]
        db_out = (beta_out[0] + beta_out[1] - beta_out[2] - beta_out[3]) / 2
    elif polarity[2] == polarity[3]:
        beta_out[0] = beta_out[0] if polarity[0] != polarity[2] else 360 - beta_out[0]
        beta_out[1] = beta_out[1] if polarity[1] != polarity[2] else 360 - beta_out[1]
        db_out = (beta_out[1] + beta_out[2] - beta_out[0] - beta_out[3]) / 2
    elif polarity[3] == polarity[0]:
        beta_out[1] = beta_out[1] if polarity[1] != polarity[0] else 360 - beta_out[1]
        beta_out[2] = beta_out[2] if polarity[2] != polarity[0] else 360 - beta_out[2]
        db_out = (beta_out[2] + beta_out[3] - beta_out[0] - beta_out[1]) / 2
    else:
        beta_out[0] = 360 - beta_out[0]
        beta_out[1] = 360 - beta_out[1]
        db_out = (beta_out[0] + beta_out[1] - beta_out[2] - beta_out[3]) / 2

    return beta_out, abs(db_out), untouched_beta


# computes:
#  - raw delta beta: raw difference in angles
# XXX - delta beta: only defined for roughly co-planar atoms (z=constant for all atoms) necessary for frame of reference
# can fix the problem with coplanar atoms by setting one corner to the origin, rotating all points so a second corner
# is also at z=0 and then rotating so the other two corners are equally close to z=0
#  - beta_in: all atoms flattened to z=0 and then the same delta beta computation takes place
#  - beta_out:
# the points should be input with their indices (1-indexed) just following adjacent angles
# i.e.
# 1-Ag  --  5-I  --  3-Bi
#  |                  |
# 6-I                7-I
#  |                  |
# 4-Bi  --  8-I  --  2-Ag
# bond=[(4, 8, 2), (2, 7, 3), (3, 5, 1), (1, 6, 4)]
# a shiftmap can also be defined, in this case:
# shiftmap = {1: '0,0,0', 2: '0,0,0', 3: '0,0,0', 4: '0,0,0', 5: '0,0,0', 6: '0,0,0', 7: '0,0,0', 8: '0,0,0'}
def new_robust_delta_beta(readfile, bonds, shiftmap=None, method="old"):
    # NOTE: no reference frame is applied for raw angles at the moment
    atom_dict = pull_data(readfile, bonds, shiftmap)

    # computes the four raw angles
    bond_coords = [[np.array(atom_dict[x]) for x in bond] for bond in bonds]
    beta = [angle(*bond) for bond in bond_coords]

    # computes four beta-in angles
    beta_in = find_db_in(bond_coords)

    # computes four beta-out angles
    beta_out, db_out, untouched_beta = find_db_out(bond_coords)
    # beta_out = untouched_beta

    db_1 = np.abs(beta[0] + beta[1] - beta[2] - beta[3]) / 2
    db_2 = np.abs(beta[0] + beta[3] - beta[2] - beta[1]) / 2
    db_in_1 = np.abs(beta_in[0] + beta_in[1] - beta_in[2] - beta_in[3]) / 2
    db_in_2 = np.abs(beta_in[0] + beta_in[3] - beta_in[2] - beta_in[1]) / 2
    db_out_1 = np.abs(beta_out[0] + beta_out[1] - beta_out[2] - beta_out[3]) / 2
    db_out_2 = np.abs(beta_out[0] + beta_out[3] - beta_out[2] - beta_out[1]) / 2
    if method == "old":
        return max(db_1, db_2)
    if method == "in":
        return max(db_in_1, db_in_2)
    if method == "out":
        return db_out
        # return max(db_out_1, db_out_2)
    if method == "new":
        if (max(db_in_1, db_in_2) + db_out) / 2 > 5 and False:
            print(readfile)
            print(beta_in)
            print(beta_out)
            print((max(db_in_1, db_in_2) + db_out) / 2)
            print()
        return max(db_in_1, db_in_2), db_out
        # return (4 * max(db_in_1, db_in_2) + db_out) / 5
        # return (max(db_in_1, db_in_2) + max(db_out_1, db_out_2)) / 2

# 16 pos
# flip 2 if: (1) all 4 up, (2) all 4 down, (3) opposites match 1 way, (4) opposites match the other way, (5)

