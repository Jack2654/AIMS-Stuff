# This file compiles basic functions for analyzing aims.out output files
# In particular, the following functions exist here:
#
# -> calc_complete(filepath):
#               determines whether the calculation in the filepath is complete by searching for "Have a nice day."
# -> find_total_energy(aims):
#               determines the total energy of the calculation in the given path
# ->
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




