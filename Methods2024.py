# This file compiles
# In particular, the following functions exist here:
#
# ->

# imports:
import os
import time
import BasicAimsOut as bao
import BasicControl as bc
import BasicGeo as bg
import numpy.polynomial.polynomial as poly
from scipy.optimize import minimize
from ase.io import read, write
import numpy as np
import math


# given the lengths of the Ag-midpoint, midpoint-Bi segments and the desired angle (in rad), returns the horizontal
# displacement of the midpoint (I)
# derivation in yellow japanese notebook
def make_angle(distA, distB, angle):
    temp = (distA + distB) ** 2
    temp = temp / (math.tan(angle) ** 2)
    temp = temp + 4 * distA * distB
    temp = math.sqrt(temp)
    temp = ((distA + distB) / (math.tan(angle))) + temp
    temp = temp / 2
    return temp


# creates delta beta in beta average structure
def induce_delta_beta(base_structure, writefile, angle1, angle2):
    lv = bg.lattice_vectors(base_structure)
    species = [x.split()[4] for x in bg.atoms(base_structure)]
    atom = bg.atoms_trimmed(base_structure)
    # atoms 7, 8, 9, 10 (indices 6, 7, 8, 9) are the in-plane iodines to be moved
    distA = atom[6][1]
    distB = atom[3][1] - distA

    disp1 = make_angle(distA, distB, angle1)
    if angle1 > math.pi:
        disp1 = -1 * make_angle(distA, distB, 2 * math.pi - angle1)

    disp2 = make_angle(distA, distB, angle2)
    if angle2 > math.pi:
        disp2 = -1 * make_angle(distA, distB, 2 * math.pi - angle2)

    atom[6][0] = disp1
    atom[7][1] = atom[3][1] - disp2
    atom[8][0] = -disp2
    atom[9][1] = atom[3][1] + disp1

    with open(writefile, "w") as f:
        for lat in lv:
            f.write(f'lattice_vector %s\n' % (" ".join([str(x) for x in lat])))
        for count, at in enumerate(atom):
            f.write(f'atom %s %s %s %s\n' % (at[0], at[1], at[2], species[count]))

