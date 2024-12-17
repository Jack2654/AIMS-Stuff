# This file compiles basic functions for analyzing aims.out output files
# In particular, the following functions exist here:
#
# -> calc_complete(filepath):
#               determines whether the calculation in the filepath is complete by searching for "Have a nice day."
# -> find_total_energy(aims):
#               determines the total energy of the calculation in the given path
# -> find_shielding(aims):
#               determines total shielding of the first element in the output file
# ->
# ->
# ->

# imports:
import os


# functions:
def calc_complete(filepath):
    if not os.path.exists(filepath):
        return False
    complete = False
    with open(filepath, "r") as f:
        for line in f:
            if "Have a nice day." in line:
                complete = True
    return complete


def find_total_energy(aims):
    with open(aims, "r") as f:
        for line in f:
            if "| Total energy of" in line:
                return line.split()[11]
    return "ERROR: no total energy found in" + aims


def all_energies(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    data = []
    for direct in all_dir:
        if os.path.isdir(base + direct):
            output = base + direct + "/aims.out"
            # print(f'%s %f' % (direct, find_shielding(output)))
            data.append(find_total_energy(output))
    return data


def analyze_band_path(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    i = 0
    for line in lines:
        if i < 2 or i > len(lines) - 3:
            temp = line.split()
            print(str(temp[1]) + " " + str(temp[2]) + " " + str(temp[3]))
            j = 0
            triggered = False
            for val in temp:
                if j > 3:
                    if j % 2 == 0:
                        if val == "0.00000" and not triggered:
                            print(temp[j+1])
                            triggered = True
                j += 1
        i += 1


def NMR_shielding_values(output):
    with open(output, "r") as f:
        found_NMR = False
        results = []
        atoms = []
        for line in f:
            if "NMR shielding tensors" in line:
                found_NMR = True
            if "Memory report" in line:
                found_NMR = False
            if found_NMR:
                if "atom" in line:
                    atoms.append(int(line.split()[1]))
                if "Total:" in line:
                    results.append(float(line.split()[1]))
    return [(atoms[x], results[x]) for x in range(len(atoms))]


# returns all shielding values found in a given folder
def all_shieldings(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    data = []
    for direct in all_dir:
        if os.path.isdir(base + direct):
            output = base + direct + "/aims.out"
            # print(f'%s %f' % (direct, find_shielding(output)))
            data.append(NMR_shielding_values(output)[0])
    return data


# creates shieldings.out file
def create_shield_out(base):
    with open(base + "aims.out", "r") as f:
        found_NMR = False
        values = []
        temp = [0, 0, 0]
        for line in f:
            if "NMR shielding tensors" in line:
                found_NMR = True
            if "Memory report" in line:
                if temp[0]:
                    values.append(temp)
                found_NMR = False
            if found_NMR:
                if "atom" in line:
                    if temp[0]:
                        values.append(temp)
                        temp = [0, 0, 0]
                    temp[0] = int(line.split()[1])
                    temp[1] = line.split()[2].replace("(", "").replace(")", "").replace(":", "")
                if "Total:" in line:
                    temp[2] = float(line.split()[1])
    with open(base + "shieldings.out", "w") as f:
        f.write("# Computed total shielding values in ppm\n")
        f.write("# Often results are referenced against TMS (~31.1846ppm)\n")
        f.write("#\t Atom\t Species\t Shielding (ppm)\n")
        for val in values:
            f.write(f'\t %s\t %s\t\t %s\n' % (val[0], val[1], val[2]))
    shieldings = [val[2] for val in values]
    print(base)
    print(" ".join([str(x) for x in shieldings]))
