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

from scipy.spatial.transform import Rotation as R
from ase.io import read, write
import VectorToolkit as vt
import numpy as np
from numpy import pi


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
                #if x > lat[i] + lv[i][i]:
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

    for line in open(filepath + "/geometry.in"):
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
        mag = math.sqrt((vec[0])**2 + (vec[1])**2 + (vec[2])**2)
        mags.append(mag)
        # print(mag)
    # alpha:
    alpha = (180 / np.pi) * np.arccos((lat[1][0] * lat[2][0] + lat[1][1] * lat[2][1] + lat[1][2] * lat[2][2]) / (mags[1] * mags[2]))
    print(alpha)
    # beta:
    beta = (180 / np.pi) * np.arccos((lat[0][0] * lat[2][0] + lat[0][1] * lat[2][1] + lat[0][2] * lat[2][2]) / (mags[0] * mags[2]))
    print(beta)
    # gamma
    gamma = (180 / np.pi) * np.arccos((lat[0][0] * lat[1][0] + lat[0][1] * lat[1][1] + lat[0][2] * lat[1][2]) / (mags[0] * mags[1]))
    print(gamma)


def delta_d(filename, center, edges, debug=False):
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
    result *= (1/6) * (1 / (avg_distance ** 2))
    if debug:
        print('atom %.5f %.5f %.5f M' % (at_center[0], at_center[1], at_center[2]))
        for atom in at_edges_moved:
            print('atom %.5f %.5f %.5f H' % (atom[0], atom[1], atom[2]))
    return result


def sigma_squared(filename, center, edges, debug=False):
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


def angle(pt1, center, pt2):
    a = np.array(pt1)
    b = np.array(center)
    c = np.array(pt2)
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    result = np.arccos(cosine_angle)
    return np.degrees(result)


def angle_info(readfile, pts, shiftmap):
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
    b = math.sqrt((connect_vect[0]**2) / (connect_vect[0]**2 + connect_vect[1]**2))
    a = math.sqrt(1 - b**2)
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
        # res.append(angle(p1, in_plane_point, p3))
        # res.append(angle(p1, out_of_plane_point, p3))
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
    b = math.sqrt((connect_vect[0]**2) / (connect_vect[0]**2 + connect_vect[1]**2))
    a = math.sqrt(1 - b**2)
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


def disturb_positions(readfile, writefile, max_disturbance):
    lv = lattice_vectors(readfile)
    at = atoms(readfile)
    with open(writefile, "w") as f:
        for lat in lv:
            f.write('lattice_vector\t%.8f\t%.8f\t%.8f\n' % (lat[0], lat[1], lat[2]))
        for atom in at:
            temp = atom.split()
            if "Cs" in temp:
                f.write('atom\t%s\t%s\t%s\t%s\n' % (temp[1], temp[2], temp[3], temp[-1]))
                continue
            coords = [float(x) for x in temp[1:4]]
            for i in range(3):
                coords[i] += max_disturbance * ((2 * random.random()) - 1)
            f.write('atom\t%.8f\t%.8f\t%.8f\t%s\n' % (coords[0], coords[1], coords[2], temp[-1]))


def maurer_displacement(readfile, corners, shiftmap, center):
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
    center_actual = atoms[center-1]
    displacement = [centroid[i] - center_actual[i] for i in range(3)]
    displacement[2] = 0
    diag_one = [1 / math.sqrt(2), 1 / math.sqrt(2), 0]
    diag_two = [1 / math.sqrt(2), -1 / math.sqrt(2), 0]

    # print(f'Centroid:\t\t%0.5f %0.5f %0.5f' % (centroid[0], centroid[1], centroid[2]))
    # print(f'Center:\t\t\t%0.5f %0.5f %0.5f' % (center_actual[0], center_actual[1], center_actual[2]))
    # print(f'Displacement:\t%0.5f %0.5f %0.5f' % (displacement[0], displacement[1], displacement[2]))

    proj_one = abs(np.dot(displacement, diag_one))
    proj_two = abs(np.dot(displacement, diag_two))

    # code to verify projections produce original result, remove abs()
    # diag_one = [x * proj_one for x in diag_one]
    # diag_two = [x * proj_two for x in diag_two]
    # total = [(diag_one[i] + diag_two[i]) * np.linalg.norm(displacement) for i in range(3)]
    # print(total)
    return proj_one, proj_two
