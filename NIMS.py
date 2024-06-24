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

name = "six"
base = f'../../Documents/NIMS/nano/Pd/raw_inputs/%s.xyz' % name
write = f'../../Documents/NIMS/nano/Pd/pure/%s.xyz' % name
# change_elem_xyz(base, write, "Pd")

name = "six"
base = f'../../Documents/NIMS/nano/Pd/pure/%s.xyz' % name
write = f'../../Documents/NIMS/nano/Pd/pure/%s.in' % name
xyz_to_cq(base, write, lattice=(39, 39, 39), move=True)

base = "../../Documents/NIMS/tmp/Pd/volume/base"
dest_start = "../../Documents/NIMS/tmp/Pd/volume/"
# values = [["111", (1, 1, 1)], ["222", (2, 2, 2)], ["333", (3, 3, 3)], ["444", (4, 4, 4)], ["555", (5, 5, 5)],
#           ["666", (6, 6, 6)], ["777", (7, 7, 7)], ["888", (8, 8, 8)]]
# values = [["025", 25], ["050", 50], ["075", 75], ["100", 100], ["125", 125], ["150", 150]]
values = ["3.75", "3.76", "3.77", "3.78", "3.79", "3.80", "3.81", "3.82", "3.83", "3.84",
          "3.85", "3.86", "3.87", "3.88", "3.89", "3.90", "3.91", "3.92", "3.93", "3.94",
          "3.95", "3.96", "3.97", "3.98", "3.99", "4.00", "4.01", "4.02", "4.03", "4.04"]
temp = 3.75
for val in values:
    dest = dest_start + val
    # make_calc(base, dest, cutoff=100, kpt=(8, 8, 8), lattice=temp)
    temp += 0.01


base = "../../Documents/NIMS/tmp/Au/DOS/DOS.dat"
with open(base, "r") as f:
    energies = []
    DOS = []
    for line in f.readlines():
        if "#" in line or "&" in line:
            continue
        else:
            temp = line.split()
            energies.append((float(temp[0]) - 4.20489) * 0.0367493)
            DOS.append(float(temp[1]))
    plt.plot(energies, DOS)
    plt.vlines(-4.20489 * 0.0367493, 0, 15, color='r')
    plt.xlim([-0.6, 0.3])
    # plt.show()
