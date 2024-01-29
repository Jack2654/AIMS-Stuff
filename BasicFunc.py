# This file compiles basic functions for doing os operations
# In particular, the following functions exist here:
#
# -> run_AIMS(directory):
#                   runs aims in the given directory
# -> run_all(base, ignore=True);
#                   runs aims in every directory in the given base. If ignore, runs regardless of
#                   whether aims has already been run in given directory
# -> mass(element):
#                   returns the nominal mass of the given element
# -> create_min_s_defaults(filepath):
#                   creates a folder of min+s defaults in the given filepath
# -> fit_poly(x, y, order):
#                   fits a polynomial of the order specified to the given data and
#                   returns the coefficients in increasing order
# -> minimize_func(coeff, guess=0.00):
#                   returns the optimized minimum of the given function specified by the coefficients
# -> merge_files(file1, file2, destination):
#                   merges two files interactively line by line
# ->

# imports:
import os
import time
import BasicAimsOut as bao
from molmass import Formula
from mins import Species
import numpy.polynomial.polynomial as poly
from scipy.optimize import minimize
from ase.io import read, write
import numpy as np


# functions:
def run_AIMS(directory):
    if not os.path.exists(directory + "geometry.in"):
        print("Bad path given to run_AIMS" + directory + "control.in")
    if not os.path.exists(directory + "control.in"):
        print("Bad path given to run_AIMS" + directory + "control.in")
    temp = os.getcwd()
    os.chdir(directory)
    mpirun = "/opt/homebrew/bin/orterun -N 8 "
    aims_exec = "/Users/jackmorgenstein/FHI-aims/Repository/builds/build_06_09_23/aims.230612.mpi.x"
    path = " > aims.out 2>&1"
    command = mpirun + aims_exec + path
    os.system(command)
    os.chdir(temp)


def run_all(base, ignore=False):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    os.system("ulimit -s hard")
    os.system("export TMPDIR=/tmp")
    os.system("export OMP_NUM_THREADS=1")
    for direct in all_dir:
        if os.path.isdir(base + direct):
            dir_name = direct
            directory = base + dir_name
            if not bao.calc_complete(directory + "/aims.out") or ignore:
                current = time.time()
                run_AIMS(directory + "/")
                print("Completed calculation for " + directory + " in " + str(time.time() - current) + " seconds")


def mass(element):
    return Formula(element).nominal_mass


def create_min_s_defaults(filepath):
    elem_dict = {'H': '1', 'He': '2', 'Li': '3', 'Be': '4', 'B': '5', 'C': '6', 'N': '7', 'O': '8',
                 'F': '9', 'Ne': '10', 'Na': '11', 'Mg': '12', 'Al': '13', 'Si': '14', 'P': '15',
                 'S': '16', 'Cl': '17', 'Ar': '18', 'K': '19', 'Ca': '20', 'Sc': '21', 'Ti': '22',
                 'V': '23', 'Cr': '24', 'Mn': '25', 'Fe': '26', 'Co': '27', 'Ni': '28', 'Cu': '29',
                 'Zn': '30', 'Ga': '31', 'Ge': '32', 'As': '33', 'Se': '34', 'Br': '35', 'Kr': '36',
                 'Rb': '37', 'Sr': '38', 'Y': '39', 'Zr': '40', 'Nb': '41', 'Mo': '42', 'Tc': '43',
                 'Ru': '44', 'Rh': '45', 'Pd': '46', 'Ag': '47', 'Cd': '48', 'In': '49', 'Sn': '50',
                 'Sb': '51', 'Te': '52', 'I': '53', 'Xe': '54', 'Cs': '55', 'Ba': '56', 'La': '57',
                 'Ce': '58', 'Pr': '59', 'Nd': '60', 'Pm': '61', 'Sm': '62', 'Eu': '63', 'Gd': '64',
                 'Tb': '65', 'Dy': '66', 'Ho': '67', 'Er': '68', 'Tm': '69', 'Yb': '70', 'Lu': '71',
                 'Hf': '72', 'Ta': '73', 'W': '74', 'Re': '75', 'Os': '76', 'Ir': '77', 'Pt': '78',
                 'Au': '79', 'Hg': '80', 'Tl': '81', 'Pb': '82', 'Bi': '83', 'Po': '84', 'At': '85',
                 'Rn': '86', 'Fr': '87', 'Ra': '88', 'Ac': '89', 'Th': '90', 'Pa': '91', 'U': '92',
                 'Np': '93', 'Pu': '94', 'Am': '95', 'Cm': '96', 'Bk': '97', 'Cf': '98', 'Es': '99',
                 'Fm': '100', 'Md': '101', 'No': '102'}
    for x in elem_dict.keys():
        species = Species(x)
        species.write_file(filepath)


def fit_poly(x, y, order):
    output = poly.polyfit(x, y, deg=order)
    return output
    # printed = "The fit is: "
    # printed += str(output[0])
    # i = 0
    # for elem in output:
    # if i != 0:
    # printed += " + " + str(elem) + " x^" + str(i)
    # i += 1
    # print(printed)


def minimize_func(coeff, guess=0.90, bnds=None):
    fun_fit = poly.Polynomial(coeff)
    return minimize(fun_fit, x0=guess, bounds=bnds).x[0]


def maximize_func(coeff):
    fun_fit = poly.Polynomial(-1*coeff)
    return minimize(fun_fit, x0=0).x[0]


def copy_to_timewarp(source, destination):      # works :D
    command = "scp -r " + source + " jhm48@timewarp.egr.duke.edu:" + destination
    os.system(command)


def mkdir_tw():
    os.system("ssh jhm48@timewarp.egr.duke.edu /bin/bash")
    os.system("mkdir brother_py")


def merge_files(file1, file2, writefile):
    with open(file1, "r") as f:
        f1 = f.readlines()
    with open(file2, "r") as f:
        f2 = f.readlines()
    max1 = len(f1)
    max2 = len(f2)
    counter1 = 0
    counter2 = 0
    with open(writefile, "w") as f:
        while True:
            if counter1 > max1 and counter2 > max2:
                break

            if counter1 < max1:
                current1 = f1[counter1]
            else:
                f.writelines(f2[counter2:len(f2)-1])
                break

            if counter2 < max2:
                current2 = f2[counter2]
            else:
                f.writelines(f1[counter1:len(f2) - 1])
                break

            # print("Comparison: " + current1.strip() + " to " + current2.strip())

            if current1 == current2:
                f.write(current1)
            else:
                processed = next_match(f1, f2, counter1, counter2)
                if processed == "same":
                    f.write(current1)

                elif processed.startswith("1"):
                    index = int(processed.split(" ")[3])
                    if index == counter2 + 1 or index - counter2 > 100:
                        f.write(current1)
                        if index == counter2 + 1:
                            counter2 = index
                    else:
                        print("File 1, line " + str(counter1 + 1) + ":")
                        print(current1)
                        print("File2, lines " + str(counter2 + 2) + " to " + str(index + 1) + ":")
                        print("".join(f2[counter2 + 1:index + 1]))
                        inp = int(input("Press 1 for file 1, 2 for file 2: "))
                        if inp == 1:
                            f.write(current1)
                            counter2 = index
                        elif inp == 2:
                            f.write("".join(f2[counter2:index + 1]))
                            counter2 = index
                        else:
                            print("ERROR")

                elif processed.startswith("2"):
                    index = int(processed.split(" ")[3])
                    if index == counter1 + 1 or index - counter1 > 100:
                        f.write(current1)
                        if index == counter1 + 1:
                            counter1 = index
                    else:
                        print("File1, lines " + str(counter1 + 2) + " to " + str(index + 1) + ":")
                        print("".join(f1[counter1 + 1:index + 1]))
                        print("File 2, line " + str(counter2 + 1) + ":")
                        print(current2)
                        inp = int(input("Press 1 for file 1, 2 for file 2: "))
                        if inp == 1:
                            f.write(current1)
                            counter1 = index
                        elif inp == 2:
                            f.write("".join(f2[counter2:index + 1]))
                            counter1 = index
                        else:
                            print("ERROR")

                elif processed == "no match":
                    print("File1:")
                    print(current1)
                    print("File2:")
                    print(current2)
                    inp = int(input("Press 1 for file 1, 2 for file 2, 3 for both: "))
                    if inp == 1:
                        f.write(current1)
                    elif inp == 2:
                        f.write(current2)
                    elif inp == 3:
                        f.write(current1)
                        f.write(current2)
                    else:
                        print("ERROR")

                elif processed == "unique 1":
                    f.write(current1)
                    counter2 -= 1

                elif processed == "unique 2":
                    f.write(current2)
                    counter1 -= 1

                else:
                    print(processed)
                    print("c1: " + str(counter1) + ", c2: " + str(counter2))
                    print("both have matches: " + current1.strip() + " versus " + current2.strip())
                    inp = int(input("1 for put 1, hold 2, 2 for put 2, hold 1: "))
                    if inp == 1:
                        f.write(current1)
                        counter2 -= 1
                    if inp == 2:
                        f.write(current2)
                        counter1 -= 1
            counter1 += 1
            counter2 += 1
    with open(writefile, "r") as f:
        print(f.readlines())


def next_match(f1, f2, c1, c2):
    line1 = f1[c1].strip()
    line2 = f2[c2].strip()
    if line1 == line2:
        return "same"

    if line1 == "":
        index = match_exists(line2, f1, c1)
        if not index == 0:
            return "2 match at " + str(index)
        return "unique 2"

    if line2 == "":
        index = match_exists(line1, f2, c2)
        if not index == 0:
            return "1 match at " + str(index)
        return "unique 1"

    index_1_match = match_exists(line1, f2, c2)
    index_2_match = match_exists(line2, f1, c1)

    if index_1_match == index_2_match == 0:
        return "no match"
    if index_1_match == 0:
        return "unique 1"
    if index_2_match == 0:
        return "unique 2"

    return "e: " + str(index_1_match) + " " + str(index_2_match)



def match_exists(base, lines, counter):
    i = counter
    while i < len(lines):
        current = lines[i].strip()
        if current == base:
            return i
        i += 1
    return 0


def disturb_geos():
    header = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/setup/"
    file_input = header + "geometry_trimmed.cif"
    displacements = ["0.05", "0.10", "0.15", "0.20"]
    for disp in displacements:
        for i in range(10):
            cif_out = header + "geo_" + disp[-2:] + "_" + str(i) + ".cif"
            command = "atomsk " + file_input + " -disturb " + disp + " " + cif_out
            os.system(command)
            with open(cif_out, "r") as f:
                lines = f.readlines()
            with open(cif_out, "w") as f:
                f.write("data_image0\n")
                f.writelines(lines)
            a = read(cif_out)
            a.write(cif_out[:-4] + ".in")


def mkdirs_and_move():
    header = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/"
    displacements = ["05", "10", "15", "20"]
    for disp in displacements:
        command = "mkdir " + header + disp + "/"
        for i in range(10):
            os.system(command + str(i))
            command2 = "mv " + header + disp + "/geo_" + disp + "_" + str(i) + ".in "
            command2 += header + disp + "/" + str(i) + "/geometry.in"
            os.system(command2)


def put_control_files():
    header = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/"
    displacements = ["05", "10", "15", "20"]
    for disp in displacements:
        for i in range(10):
            command = "cp " + header + "setup/control.in "
            command += header + disp + "/" + str(i) + "/."
            os.system(command)


def read_BZ_corners(filename):
    with open(filename, "r") as f:
        lns = f.readlines()
    ret = []
    for ln in lns:
        temp = ln.split()
        temp = [x.replace("(", "") for x in temp]
        temp = [x.replace(")", "") for x in temp]
        temp = [x.replace(",", "") for x in temp]
        ret.append([float(x) for x in temp[0:3]])
    return ret
