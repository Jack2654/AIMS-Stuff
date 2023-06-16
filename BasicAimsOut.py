# This file compiles basic functions for analyzing aims.out output files
# In particular, the following functions exist here:
#
# -> run_AIMS(directory):
#                   runs aims in the given directory
# -> run_all(base, ignore=True);
#                   runs aims in every directory in the given base. If ignore, runs regardless of
#                   whether aims has already been run in given directory
# ->
# ->
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
