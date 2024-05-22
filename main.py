import BasicFunc
import Geometry
import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt
import BasicBandOut as bbo
import BasicAimsOut as bao
import BasicControl as bc
import BasicFunc as bf
import matplotlib.pyplot as plt
import math
import random
import time
import numpy as np
import os


settings = "../../FHI-aims/Double_Perovskites/New_Structures/random/bands_look_good/180_old_model/180_new2/settings_final.in"
# pt.mulliken_plot(settings, debug=False, quiet=False, save=True)

# atom dict for DMF
atom_dict = [[9, 10, 11], [6, 7, 8], [12]]
color_dict = ['b', 'r', 'g']
# pt.NMR_average("../../FHI-aims/French_NMR/DMF/Shield_MD/", atom_dict, color_dict)

NH2 = [1, 116, 2, 105, 3, 96, 4, 85, 5, 78, 6, 71]
CH2 = [63, 64, 103, 104, 97, 98, 83, 84, 79, 80, 69, 70]
Benzene = [65, 66, 67, 68, 99, 100, 101, 102, 61, 62, 81, 82]
Cyclohexane = [106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                    72, 73, 74, 75, 76, 77, 117, 118, 119, 120]
atom_dict = [NH2, CH2, Benzene, Cyclohexane]
# atom dict for TA
atom_dict[0] = [116, 105, 96, 85, 78, 71]
atom_dict = [[x - 12 for x in part] for part in atom_dict]
color_dict = ['r'] * 6 + \
             ['c'] * 6 + \
             ['b'] * 3 + \
             ['g'] * 3
new_atom_dict = [[atom_dict[0][0]], [atom_dict[0][1]], [atom_dict[0][2]],
                 [atom_dict[0][3]], [atom_dict[0][4]], [atom_dict[0][5]],

                 atom_dict[1][0:2], atom_dict[1][2:4], atom_dict[1][4:6],
                 atom_dict[1][6:8], atom_dict[1][8:10], atom_dict[1][10:],

                 atom_dict[2][0:4], atom_dict[2][4:8], atom_dict[2][8:],

                 atom_dict[3][0:10], atom_dict[3][10:20], atom_dict[3][20:]]


# 1 NH group, 6 CH2 groups, 3 Benzene groups, 3 Cyclohexane groups
folder = "../../FHI-aims/French_NMR/NMR_results/TA_NMR/"
# pt.NMR_density(folder, new_atom_dict, color_dict, average=True, width=0.01)
folder = "../../FHI-aims/French_NMR/TA_MD/MD_NMR/Shield_00/"
pt.NMR_density(folder, new_atom_dict, color_dict, average=True, width=0.01)

