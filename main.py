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
directories = ["20.0_hse/"]
for count, angle in enumerate(directories):
    folder = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/" + str(base_angle) + "/" + directories[count]
    # print(folder)
    res = []
    res.append(bbo.band_info(folder, "band1001.out", spin_splitting=True, verbose=False))
    res.append(bbo.band_info(folder, "band1002.out", spin_splitting=True, verbose=False))
    res.append(bbo.band_info(folder, "band1003.out", spin_splitting=True, verbose=False))
    res.append(bbo.band_info(folder, "band1004.out", spin_splitting=True, verbose=False))
    # res.append(bbo.band_info(folder, "band1005.out", spin_splitting=True, verbose=False))
    print(" ".join(res))
    pt.mulliken_plot(folder + "settings.in", debug=False, quiet=False, save=False)

read_folder = "../../FHI-aims/French_NMR/random/TMS_work/TMS_MD/MD/geometries/geometry"
write_folder = "../../FHI-aims/French_NMR/random/TMS_work/TMS_MD/NMR_dense/"
# bf.make_MD_DMF(read_folder, write_folder, calcs=10, step=10, species=["C", "H"])

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
# structures to fix:
# AgBi:

# Pb:

# TlBi:


file = "../../FHI-aims/Double_Perovskites/New_Structures/Pb/bands/data.txt"
# pt.plot_beta_avg_corr(file)
file = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/data.txt"
# pt.plot_beta_avg_corr(file)

read = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/height/180_diff/geometry.in"
write = read
# ots.constrain_base_relax(read, write)


