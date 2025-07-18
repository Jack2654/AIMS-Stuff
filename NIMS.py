import os
import matplotlib.pyplot as plt
import numpy as np


def make_calc(base, dest, cutoff=100, kpt=(1, 1, 1), lattice=None):
    command = f'rm -r %s' % dest
    os.system(command)
    command = f'cp -r %s %s' % (base, dest)
    os.system(command)
    with open(f'%s/Conquest_input' % dest, "r") as f:
        inp = f.readlines()
    with open(f'%s/Conquest_input' % dest, "w") as f:
        for line in inp:
            temp = line
            if "Grid.GridCutoff" in line:
                temp = f'Grid.GridCutoff           %s\n' % cutoff
            elif "Diag.MPMeshX" in line:
                temp = f'Diag.MPMeshX              %s\n' % kpt[0]
            elif "Diag.MPMeshY" in line:
                temp = f'Diag.MPMeshY              %s\n' % kpt[1]
            elif "Diag.MPMeshZ" in line:
                temp = f'Diag.MPMeshZ              %s\n' % kpt[2]
            f.write(temp)
    if lattice:
        with open(f'%s/geometry.in' % dest, "r") as f:
            inp = f.readlines()
        with open(f'%s/geometry.in' % dest, "w") as f:
            f.write(f'%s   0   0\n' % lattice)
            f.write(f'0   %s   0\n' % lattice)
            f.write(f'0   0   %s\n' % lattice)
            for line in inp:
                if len(line.split()) == 3:
                    continue
                else:
                    f.write(line)


def change_elem_xyz(base, write, elem):
    with open(base, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        for line in lines:
            temp = line.split()
            if len(temp) == 4:
                f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (elem, temp[1], temp[2], temp[3]))
            else:
                f.write(line)


def xyz_to_cq(base, write, lattice=(0, 0, 0), move=False):
    movem = "F F F"
    if move:
        movem = "T T T"
    with open(base, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        f.write(f'%s 0 0\n' % lattice[0])
        f.write(f'0 %s 0\n' % lattice[1])
        f.write(f'0 0 %s\n' % lattice[2])
        elem_dict = {}
        elem_count = 1
        for line in lines:
            temp = line.split()
            if len(temp) == 1:
                f.write(line.strip() + "\n")
            if len(temp) == 4:
                if temp[0] not in elem_dict:
                    elem_dict[temp[0]] = elem_count
                    elem_count += 1
                f.write(f'%s\t\t%s %s\n' % (" ".join(temp[1:]), elem_dict[temp[0]], movem))


def cq_to_aims(base, write, frac=False):
    extra = "_frac" if frac else ""
    with open(base, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        f.write(f'lattice_vector %s' % lines[0])
        f.write(f'lattice_vector %s' % lines[1])
        f.write(f'lattice_vector %s' % lines[2])
        # for line in lines[4:-2]:
        for line in lines[4:]:
            temp = line.split()
            if temp[3] == "1" and "Pd" in base:
                f.write(f'atom{extra} %s\t\t Pd\n' % " ".join(temp[:3]))
            else:
                f.write(f'atom{extra} %s\t\t Ag\n' % " ".join(temp[:3]))
            # f.write(f'atom{extra} %s\t\t Ag\n' % " ".join(temp[:3]))
        # temp = lines[-2].split()
        # f.write(f'atom{extra} %s\t\t H\n' % " ".join(temp[:3]))
        # temp = lines[-1].split()
        # f.write(f'atom{extra} %s\t\t H\n' % " ".join(temp[:3]))

    # print(f'Transformed Input File: {base}')
    # print(f'Into aims file: {write}')


def bohr_to_a(base, frac=False, convert=True):
    extra = "_frac" if frac and not convert else ""
    with open(base, "r") as f:
        lines = f.readlines()
    with open(base, "w") as f:
        for line in lines:
            temp = line.split()
            if "lattice" in temp[0]:
                val = [0.529177249 * float(x) for x in temp[1:]]
                lat = max(val)
                f.write(f'lattice_vector %s %s %s\n' % (val[0], val[1], val[2]))
            elif "atom" in temp[0]:
                val = [0.529177249 * float(x) for x in temp[1:4]] if not frac else [float(x) for x in temp[1:4]]
                if convert:
                    val[0] *= lat
                    val[1] *= lat
                    val[2] *= lat
                # for elem in range(len(val)):
                #     if val[elem] > lat / 2:
                #         val[elem] = val[elem] - lat
                #     val[elem] += lat/2
                f.write(f'atom{extra} %s %s %s %s\n' % (val[0], val[1], val[2], temp[4]))
            else:
                print("Unexpected:")
                print(line)
    # print(f'Translated Bohr to Angstrom: {base}')
    # print()


def a_to_bohr(base, write):
    with open(base, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        counter = 0
        for line in lines:
            temp = line.split()
            if counter == 3:
                f.write(line)
            else:
                vals = [float(x) * 1.8897259886 for x in temp[0:3]]
                vals = [str(x) for x in vals[0:3]]
                f.write(f'%s %s\n' % (" ".join(vals), " ".join(temp[3:])))
                # 1.8897259886
            counter += 1


def distance(atom1, atom2):
    return np.linalg.norm(np.array(atom1) - np.array(atom2))


def nn_dist(base):
    with open(base, "r") as f:
        lines = f.readlines()
    atoms = [[float(x) for x in line.split()[1:4]] for line in lines if "atom" in line]
    np_atoms = atoms[:923]
    ad_atoms = atoms[923:]
    dist = 100
    for ad in ad_atoms:
        for atom in np_atoms:
            cur_dist = distance(ad, atom)
            if cur_dist < dist:
                dist = cur_dist
    return dist


def gen_nn_dists():
    base_path = "../../Documents/NIMS/Adsorbate_Results/Info_06_03/coord_next/"
    read_folder = os.path.join(base_path, "Ag6")
    files = [file for file in os.listdir(read_folder) if "opt" in file and "aims" not in file and "xyz" not in file]
    for adsorbate in ["H2", "OH", "_H_", "_O_"]:
        base_files = [file for file in files if adsorbate in file and "fixed" not in file]
        fixed_files = [file for file in files if adsorbate in file and "fixed" in file]
        dest_base = os.path.join(base_path, "plotting/Ag6_base/")
        # copy over all OG geometry files
        for geo_file in base_files:
            os.system(f'cp {os.path.join(read_folder, geo_file)} {dest_base + geo_file}')
        # copy over 'fixed' geometry files, replacing non-fixed versions
        for geo_file in fixed_files:
            os.system(f'cp {os.path.join(read_folder, geo_file)} {dest_base + geo_file.replace("fixed_", "")}')

    read_path = os.path.join(base_path, "plotting/Ag6_base")
    write_path = os.path.join(base_path, "plotting/Ag6_aims")
    coords = os.listdir(read_path)
    for adsorbate in ["H2", "OH", "_H_", "_O_"]:
        cur_coords = [file for file in coords if adsorbate in file]
        cur_coords.sort()
        for file in cur_coords:
            cq_to_aims(os.path.join(read_path, file), os.path.join(write_path, "aims_" + file), frac=True)
            bohr_to_a(os.path.join(write_path, "aims_" + file), frac=True, convert=True)
            dist = nn_dist(os.path.join(write_path, "aims_" + file))
            print(f'{file.replace("coord_AgNP_", "").replace("_opt.in", "")}: {dist}')


def adsorbate_plot(data_file):
    with open(data_file, "r") as f:
        lines = f.readlines()
    colors = {'H': 'r', 'OH': 'lightblue', 'H2': 'g', 'O': 'y'}
    shapes = {'bridge': 'o', 'top': 's', 'hollow': 'D'}
    for line in lines:
        cur = line.split()
        attributes = cur[0].split("_")
        cur_dist = float(cur[1])
        cur_energy = float(cur[2])
        plt.scatter(cur_dist, cur_energy, s=25, facecolors=colors[attributes[0]], marker=shapes[attributes[2]],
                    edgecolors='black', linewidths=1.5)

    plt.plot(0, 0, color='r', label='H')
    plt.plot(0, 0, color='lightblue', label='OH')
    plt.plot(0, 0, color='g', label='H2')
    plt.plot(0, 0, color='y', label='O')
    plt.scatter(0, 0, marker='o', color='k', label='Bridge')
    plt.scatter(0, 0, marker='D', color='k', label='Hollow')
    plt.scatter(0, 0, marker='s', color='k', label='On-top')
    plt.xlim([1.5, 3.7])
    plt.xlabel("Nearest Neighbor Distance (Å)")
    plt.ylabel("Adsorption Energy (Ha)")
    plt.title("Adsorbates on Ag6 NP")
    plt.legend()
    plt.tight_layout()
    plt.show()


adsorbate_plot("../../Documents/NIMS/Adsorbate_Results/Info_06_03/coord_next/plotting/data.in")

path = "/Users/jackmorgenstein/Documents/NIMS/new_calcs/geo_files/"
for filename in os.listdir(path):
    # for filename in ["6L"]:
    continue
    folder = os.path.join(path, filename)
    for file in os.listdir(folder):
        write_path = os.path.join(folder, "aims_") + file + ".in"
        cq_to_aims(os.path.join(folder, file), write_path, frac="False")
        bohr_to_a(write_path, frac=False)

cur_path = f'{path}coord_AgNP_H2_100_top_c7_opt.in'
write = f'{path}aims_coord_AgNP_H2_100_top_c7_opt.in'
# cq_to_aims(cur_path, write, frac="True")
# bohr_to_a(write, frac="True")
# faces = ["100", "110", "111"]
# locations = "hollow", "bridge", "top"

# order of hollow, bridge, top for 100, 110, 111 for O adsorbate
particle = "O"
label_dict = [[[i for i in range(1, 7)], [i for i in range(1, 10)], [1, 5, 6, 7, 8, 9, 10]],
              [[], [1, 3], [1, 2, 3]],  # 110_bridge_b2 is missing from input files but should exist
              [[i for i in range(1, 10)], [1, 2, 3, 4, 6, 7, 8, 9], [5, 6, 7]]]  # 111_bridge_b5 should exist

# order of hollow, bridge, top for 100, 110, 111 for H2 adsorbate
particle = "H2"
label_dict = [[[i for i in range(1, 7)], [i for i in range(1, 10)], [1, 5, 6, 7, 8, 9, 10]],
              [[], [1, 2, 3], [1, 2, 3]],
              [[i for i in range(1, 10)], [i for i in range(1, 10)], [5, 6, 7]]]

# for i, set_o_labels in enumerate(label_dict):
#     base_path = f'{path}coord_AgNP_{particle}_{faces[i]}_'
#     for j, lbls in enumerate(set_o_labels):
#         for index in lbls:
#             cur_path = f'{base_path}{locations[j]}_c{str(index)}_opt.in'
#             write = f'{path}aims_{faces[i]}_{locations[j]}_c{str(index)}.in'
#             cq_to_aims(cur_path, write, frac="True")
#             bohr_to_a(write, frac="True")


type_dict = {"hollow": [6, 0, 9], "bridge": [9, 3, 9], "top": [10, ]}

files = ["110_top_b1",
         "110_bridge_b1",
         "100_top_b10",
         "100_hollow_b3",
         "100_hollow_b5",
         "100_hollow_b6"]
for file in files:
    continue
    base = f'{path}coord_AgNP_O_{file}_opt.in'
    write = f'{path}aims_{file}.in'
    cq_to_aims(base, write, frac="True")
    bohr_to_a(write, frac="True")

# site 6
energies = [-0.1886, -0.1876, -0.1846, -0.1849, -0.1830, -0.1800]
# site 7
energies = [-0.1892, -0.1881, -0.1832, -0.1845, -0.1840, -0.1813]
# plt.scatter([0, 1, 2, 3, 4, 5], energies)
# plt.xlabel("Palladium Layers")
# plt.ylabel("Adsorption Energy (Ha)")
# plt.title("Adsorption Energy of Bridge Site 7")
# plt.tight_layout()
# plt.show()

# xyz_to_cq(base, base, lattice=(40, 40, 40), move=True)
# a_to_bohr(base, base)

# directories = ["Pd1Ag5", "Pd2Ag4", "Pd3Ag3", "Pd4Ag2", "Pd5Ag1"]
directories = ["purePd", "pureAg"]
directory = "../../Documents/NIMS/Scripts/"

layer_matrix = [[1, 0, 0, 0, 0],
                [0, 12, 0, 0, 0],
                [0, 12, 24, 0, 6],
                [0, 12, 48, 8, 24],
                [0, 12, 72, 24, 54],
                [0, 12, 96, 48, 96],
                [0, 12, 120, 80, 150]]

# layers = 2
# with open(directory + f'labels_%s.txt' % layers, "w") as f:
#     f.write("iatom\t\tNPtype\t\tposition\tilayer\n")
#     cur_layer = 0
#     atom_no = 1
#     while cur_layer <= layers:
#         for val in range(layer_matrix[cur_layer][0]):
#             f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 0, cur_layer))
#             atom_no += 1
#         for val in range(layer_matrix[cur_layer][1]):
#             f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 1, cur_layer))
#             atom_no += 1
#         for val in range(layer_matrix[cur_layer][2]):
#             f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 2, cur_layer))
#             atom_no += 1
#         for val in range(layer_matrix[cur_layer][4]):
#             f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 4, cur_layer))
#             atom_no += 1
#         for val in range(layer_matrix[cur_layer][3]):
#             f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 3, cur_layer))
#             atom_no += 1
#         cur_layer += 1


for dir in directories:
    continue
    base = "../../Documents/NIMS/nano/6L/Ag/" + dir + ".xyz"
    write = "../../Documents/NIMS/nano/6L/Ag/" + dir + ".in"
    xyz_to_cq(base, write, lattice=(40, 40, 40), move=True)
    # cq_to_aims(write, write + "aims")
    a_to_bohr(write, write)

# base = "../../Documents/NIMS/tmp/Pd/volume/base"
# dest_start = "../../Documents/NIMS/tmp/Pd/volume/"
# values = [["111", (1, 1, 1)], ["222", (2, 2, 2)], ["333", (3, 3, 3)], ["444", (4, 4, 4)], ["555", (5, 5, 5)],
#           ["666", (6, 6, 6)], ["777", (7, 7, 7)], ["888", (8, 8, 8)]]
# values = [["025", 25], ["050", 50], ["075", 75], ["100", 100], ["125", 125], ["150", 150]]
# values = ["3.75", "3.76", "3.77", "3.78", "3.79", "3.80", "3.81", "3.82", "3.83", "3.84",
#           "3.85", "3.86", "3.87", "3.88", "3.89", "3.90", "3.91", "3.92", "3.93", "3.94",
#           "3.95", "3.96", "3.97", "3.98", "3.99", "4.00", "4.01", "4.02", "4.03", "4.04"]
# temp = 3.75
# for val in values:
#     dest = dest_start + val
#     # make_calc(base, dest, cutoff=100, kpt=(8, 8, 8), lattice=temp)
#     temp += 0.01


# base = "../../Documents/NIMS/nano/Pd/DOS/DOS.dat"
# with open(base, "r") as f:
#     energies = []
#     DOS = []
#     for line in f.readlines():
#         if "#" in line or "&" in line:
#             continue
#         else:
#             temp = line.split()
#             energies.append((float(temp[0]) - 4.20489) * 0.0367493)
#             DOS.append(float(temp[1]))
#     plt.plot(energies, DOS)
#     plt.vlines(-4.20489 * 0.0367493, 0, 15, color='r')
#     plt.xlim([-0.6, 0.3])
#     # plt.show()


# with open(base, "r") as f:
#     energies = []
#     DOS = []
#     for line in f.readlines():
#         if "#" in line or "&" in line:
#             continue
#         else:
#             temp = line.split()
#             energies.append(float(temp[0]))
#             DOS.append(float(temp[1]))
#     plt.plot(energies, DOS)
#     plt.vlines(0, 0, 15, color='r')
#     plt.xlim([-7, 2])
#     plt.yticks([])
#     # plt.show()
