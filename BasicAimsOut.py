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
