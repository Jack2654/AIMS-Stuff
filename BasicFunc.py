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
# ->
# ->

# imports:
import os
import time
import BasicAimsOut as bao
from molmass import Formula
from mins import Species
import numpy.polynomial.polynomial as poly
from scipy.optimize import minimize


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
                print("Completed calculation for " + dir_name + " in " + str(time.time() - current) + " seconds")


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


def minimize_func(coeff, guess=0.90):
    fun_fit = poly.Polynomial(coeff)
    return minimize(fun_fit, x0=guess).x[0]
