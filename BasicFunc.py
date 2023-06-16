# This file compiles basic functions for doing os operations
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
import time
import BasicAimsOut as bao


# functions:
def run_AIMS(directory):
    if not os.path.exists(directory + "geometry.in"):
        print("Bad path given to run_AIMS" + directory + "control.in")
    if not os.path.exists(directory + "control.in"):
        print("Bad path given to run_AIMS" + directory + "control.in")
    temp = os.getcwd()
    os.chdir(directory)
    os.system("ulimit -s hard")
    os.system("export TMPDIR=/tmp")
    os.system("export OMP_NUM_THREADS=1")
    mpirun = "/opt/homebrew/bin/orterun -N 8 "
    aims_exec = "/Users/jackmorgenstein/FHI-aims/Repository/builds/build_06_09_23/aims.230612.mpi.x"
    path = " > aims.out 2>&1"
    command = mpirun + aims_exec + path
    os.system(command)
    os.chdir(temp)


def run_all(base, ignore=False):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for direct in all_dir:
        if os.path.isdir(base + direct):
            dir_name = direct
            directory = base + dir_name
            if not bao.calc_complete(directory + "/aims.out") or ignore:
                current = time.time()
                run_AIMS(directory + "/")
                print("Completed calculation for " + dir_name + " in " + str(time.time() - current) + " seconds")
