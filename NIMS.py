import os
import matplotlib.pyplot as plt

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


def cq_to_aims(base, write):
    with open(base, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        f.write(f'lattice_vector %s' % lines[0])
        f.write(f'lattice_vector %s' % lines[1])
        f.write(f'lattice_vector %s' % lines[2])
        for line in lines[4:]:
            temp = line.split()
            if temp[3] == "1":
                f.write(f'atom %s\t\t Pd\n' % " ".join(temp[:3]))
            else:
                f.write(f'atom %s\t\t Ag\n' % " ".join(temp[:3]))


def bohr_to_a(base):
    with open(base, "r") as f:
        lines = f.readlines()
    with open(base, "w") as f:
        for line in lines:
            temp = line.split()
            if "lattice" in temp[0]:
                val = [0.529177 * float(x) for x in temp[1:]]
                f.write(f'lattice_vector %s %s %s\n' % (val[0], val[1], val[2]))
            elif "atom" in temp[0]:
                val = [0.529177 * float(x) for x in temp[1:4]]
                f.write(f'atom %s %s %s %s\n' % (val[0], val[1], val[2], temp[4]))
            else:
                print("Unexpected:")
                print(line)


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


base = f'../../Documents/NIMS/nano/2L/starting/Ag.in'
xyz_to_cq(base, base, lattice=(40, 40, 40), move=True)
a_to_bohr(base, base)
# cq_to_aims(base, write)
# bohr_to_a(write)

# directories = ["Pd1Ag5", "Pd2Ag4", "Pd3Ag3", "Pd4Ag2", "Pd5Ag1"]
directories = ["purePd", "pureAg"]
directory = "../../Documents/NIMS/Scripts/"

layer_matrix = [[1, 0,  0,   0,  0],
                [0, 12, 0,   0,  0],
                [0, 12, 24,  0,  6],
                [0, 12, 48,  8,  24],
                [0, 12, 72,  24, 54],
                [0, 12, 96,  48, 96],
                [0, 12, 120, 80, 150]]

layers = 2
with open(directory + f'labels_%s.txt' % layers, "w") as f:
    f.write("iatom\t\tNPtype\t\tposition\tilayer\n")
    cur_layer = 0
    atom_no = 1
    while cur_layer <= layers:
        for val in range(layer_matrix[cur_layer][0]):
            f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 0, cur_layer))
            atom_no += 1
        for val in range(layer_matrix[cur_layer][1]):
            f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 1, cur_layer))
            atom_no += 1
        for val in range(layer_matrix[cur_layer][2]):
            f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 2, cur_layer))
            atom_no += 1
        for val in range(layer_matrix[cur_layer][4]):
            f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 4, cur_layer))
            atom_no += 1
        for val in range(layer_matrix[cur_layer][3]):
            f.write(f'%s\t\t%s\t\t%s\t\t%s\n' % (atom_no, 3, 3, cur_layer))
            atom_no += 1
        cur_layer += 1


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
