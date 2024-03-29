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


def find_shielding(aims):
    value = "ERROR"
    with open(aims, "r") as f:
        start = False
        for line in f:
            if "WELCOME TO MAGNETIC RESPONSE CALCULATIONS" in line:
                start = True
            if start:
                if "Total" in line:
                    value = float(line.split()[1])
                    start = False
    return value


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
        for line in f:
            if "NMR shielding tensors" in line:
                found_NMR = True
            if "Memory report" in line:
                found_NMR = False
            if found_NMR and "Total:" in line:
                results.append(float(line.split()[1]))
    return "\t\t".join([str(x) for x in results])


# returns all shielding values found in a given folder
def all_shieldings(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    data = []
    for direct in all_dir:
        if os.path.isdir(base + direct):
            output = base + direct + "/aims.out"
            # print(f'%s %f' % (direct, find_shielding(output)))
            data.append(find_shielding(output))
    return data
