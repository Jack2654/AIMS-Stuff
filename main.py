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
import Methods2024 as M24

options = ["180/", "170/", "160/", "150/", "150_low_Cs/"]
for opt in options:
    folder = "../../FHI-aims/Double_Perovskites/New_Structures/base_structures/bands/" + opt
    # print(opt)
    # bbo.band_info(folder, "band1001.out", band_gap=True)
    # pt.mulliken_plot(folder + "settings.in", debug=False, quiet=False, save=False)

base_angle = 180
directories = ["00.0/", "02.5/", "05.0/", "07.5/", "10.0/", "12.5/", "15.0/", "17.5/", "20.0/"]
directories = ["20.0_new/"]
directories = ["n_2_4_experimental/", "n_2_4_theoretical/"]
bands = ['band0.out', 'band5.out', 'band6.out', 'band1.out', 'band4.out']
for count, angle in enumerate(directories):
    continue
    # folder = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/" + str(base_angle) + "/" + directories[count]
    # folder = "../../FHI-aims/Double_Perovskites/New_Structures/random/bands_look_good/20.0_new/"
    # print(folder)
    folder = f'../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/%s' % angle
    for band in bands:
        print(bbo.band_info(folder, band, verbose=False))
        bbo.band_gap(folder + band)
    pt.mulliken_plot(folder + "settings_new_plane.in", debug=False, quiet=False, save=False)
    print()

read_folder = "../../FHI-aims/French_NMR/I3/MD/geometries/geometry"
write_folder = "../../FHI-aims/French_NMR/I3/NMR/"
# bf.make_MD_DMF(read_folder, write_folder, calcs=10, step=200, species=["H"])

base_angle = 180
# angles = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20]
# directories = ["00.0/", "02.5/", "05.0/", "07.5/", "10.0/", "12.5/", "15.0/", "17.5/", "20.0/"]
angles = [0]
directories = ["00.0/"]
base = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/height/" + str(base_angle) + "/geometry.in.next_step"
for count, angle in enumerate(angles):
    continue
    ang1 = (base_angle + angle / 2) * math.pi / 180
    ang2 = (base_angle - angle / 2) * math.pi / 180
    write = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/" + str(base_angle) + "/" + \
            directories[count] + "geometry.in"
    print(write)
    # M24.induce_delta_beta(base, write, ang1, ang2)

    test = write
    print(round(bg.angle_info(test, (1, 7, 4), -1)[0], 2))
    print(round(bg.angle_info(test, (4, 8, 1), [[0, 0, 0], [0, 0, 0], [0, 1, 0]])[0], 2))
    print(round(bg.angle_info(test, (4, 9, 1), [[0, 0, 0], [0, 0, 0], [1, 1, 0]])[0], 2))
    print(round(bg.angle_info(test, (4, 10, 1), [[0, 0, 0], [0, 0, 0], [1, 0, 0]])[0], 2))


file = "../../FHI-aims/Double_Perovskites/New_Structures/Pb/bands/data.txt"
# pt.plot_beta_avg_corr(file)
file = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/data.txt"
# pt.plot_beta_avg_corr(file)


# file = "../../FHI-aims/Double_Perovskites/New_Structures/random/bands_look_good/20.0_new/band1001.out"
file = f'../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/n_4_experimental/'
bands = ["band0.out", "band5.out", "band6.out", "band1.out", "band3.out"]
for band in bands:
    res = bbo.band_gap(file, band, valence_offset=0, conduction_offset=0, display=True)[0:2]
    print(" ".join([str(x) for x in res]))
    res = bbo.band_gap(file, band, valence_offset=0, conduction_offset=2, display=True)[0:2]
    # res = bbo.band_gap(file, band, valence_offset=2, conduction_offset=0, display=True)[0:2]
    print(" ".join([str(x) for x in res]))
    # print(bbo.band_info(file, band, effective_mass=True, verbose=False))
    # print()


#################### Info for TA + 3H + 3I3 ####################
# 1H info
folder = "../../FHI-aims/French_NMR/step3/I3/NMR/"
NH2 = [10, 116, 11, 96, 12, 78, 71, 85, 105]
CH2 = [63, 64, 103, 104, 97, 98, 83, 84, 79, 80, 69, 70]
Benzene = [65, 66, 67, 68, 99, 100, 101, 102, 61, 62, 81, 82]
Cyclohexane = [106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
               72, 73, 74, 75, 76, 77, 117, 118, 119, 120]
atom_dict = [NH2, CH2, Benzene, Cyclohexane]
atom_dict = [[atom_dict[0][0:2], atom_dict[0][2:4], atom_dict[0][4:6],
              [atom_dict[0][6]], [atom_dict[0][7]], [atom_dict[0][8]]],

             [atom_dict[1][0:2], atom_dict[1][2:4], atom_dict[1][4:6],
              atom_dict[1][6:8], atom_dict[1][8:10], atom_dict[1][10:]],

             [atom_dict[2][0:4], atom_dict[2][4:8], atom_dict[2][8:]],

             [atom_dict[3][0:10], atom_dict[3][10:20], atom_dict[3][20:]]]
color_dict = ['r', 'c', 'b', 'g']
type = 'H'
xlim = [-10, 0]

# pt.NMR_average(folder, atom_dict, color_dict, width=0.005, xlim=xlim, type=type)

#################### Info for TA + 3H + 1I3 + 2I ####################
# 1H info
folder = "../../FHI-aims/French_NMR/step3/I3/2I_I3/NMR/"
NH2 = [10, 116, 11, 96, 12, 78, 71, 85, 105]
CH2 = [63, 64, 103, 104, 97, 98, 83, 84, 79, 80, 69, 70]
Benzene = [65, 66, 67, 68, 99, 100, 101, 102, 61, 62, 81, 82]
Cyclohexane = [106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
               72, 73, 74, 75, 76, 77, 117, 118, 119, 120]
atom_dict = [NH2, CH2, Benzene, Cyclohexane]
atom_dict = [[val - 4 for val in sset] for sset in atom_dict]
atom_dict = [[atom_dict[0][0:2], atom_dict[0][2:4], atom_dict[0][4:6],
              [atom_dict[0][6]], [atom_dict[0][7]], [atom_dict[0][8]]],

             [atom_dict[1][0:2], atom_dict[1][2:4], atom_dict[1][4:6],
              atom_dict[1][6:8], atom_dict[1][8:10], atom_dict[1][10:]],

             [atom_dict[2][0:4], atom_dict[2][4:8], atom_dict[2][8:]],

             [atom_dict[3][0:10], atom_dict[3][10:20], atom_dict[3][20:]]]
color_dict = ['r', 'c', 'b', 'g']
type = 'H'
xlim = [-14, 0]

# pt.NMR_average(folder, atom_dict, color_dict, width=0.005, xlim=xlim, type=type)
