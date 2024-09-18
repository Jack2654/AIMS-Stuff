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

file = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/n_2_4_experimental/"
alt_geo = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/plane_geos/n_2_4_theoretical/geometry.in"
# pt.mulliken_plot(file + "settings_new_plane.in", alt_geo=alt_geo, debug=False, save=False)

read_folder = "../../FHI-aims/French_NMR/new_multi/I3/MD/mix/geometries/geometry"
write_folder = "../../FHI-aims/French_NMR/new_multi/NMR/I3/mix/"
# bf.make_MD_DMF(read_folder, write_folder, calcs=10, step=1000, species=["H"])

file = "../../FHI-aims/Double_Perovskites/New_Structures/Pb/bands/data.txt"
# pt.plot_beta_avg_corr(file)
file = "../../FHI-aims/Double_Perovskites/New_Structures/TlBi/bands/data.txt"
# pt.plot_beta_avg_corr(file)

file = f'../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/n_4_experimental/'
file = f'../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/n_4_theoretical/'
# n_2_4 bands = ["band0.out", "band5.out", "band6.out", "band1.out", "band4.out"]
# n_3 bands = ["band0.out", "band6.out", "band7.out", "band4.out", "band3.out"]
bands = ["band0.out", "band5.out", "band6.out", "band1.out", "band3.out"]
# n_5 bands = ["band0.out", "band6.out", "band7.out", "band1.out", "band3.out"]
# n_4_6 bands = ["band1001.out", "band1004.out", "band1005.out", "band1003.out", "band1002.out"]
steps = 5
for band in bands:
    continue
    # or x in range(steps):
    res = bbo.band_gap(file, band, valence_offset=0, conduction_offset=2, steps=steps, display=True)[0:2]
    print(" ".join([str(x) for x in res]))
    # res = bbo.band_gap(file, band, valence_offset=0, conduction_offset=2, steps=steps, display=True)[0:2]
    # print(" ".join([str(x) for x in res]))

base_angle = 180
directories = ["20.0/", "02.5/", "05.0/", "07.5/", "10.0/", "12.5/", "15.0/", "17.5/", "20.0/"]
bands = ['band1001.out', 'band1002.out', 'band1003.out', 'band1004.out', 'band1005.out']
for count, angle in enumerate(directories):
    continue
    folder = "../../FHI-aims/Double_Perovskites/New_Structures/Pb/bands/" + str(base_angle) + "/" + directories[count]
    # print(folder)
    res = []
    for band in bands:
        res.append(bbo.band_info(folder, band, band_gap=False, spin_splitting=True, verbose=False))
        # print(bbo.band_gap(folder, band, valence_offset=-1, conduction_offset=1, display=True)[0])
    # pt.mulliken_plot(folder + "settings.in", debug=False, quiet=True, save=False)
    print(" ".join([str(x) for x in res]))

# pt.plot_beta_avg_corr("../../FHI-aims/Double_Perovskites/New_Structures/Pb/bands/data.txt")

# ots.ME431_lab1()

# BZ plotting code for m=4, n=4
corner_file = "../../FHI-aims/Yi/Yi_1_5_D/Results/New_Results/n_4/experimental/BZ_corners.in"
corners = np.asarray(bf.read_BZ_corners(corner_file))
adjacency = [[1, 6, 8], [7, 9], [3, 4, 10], [5, 11], [6, 5], [7], [7],
             [], [9, 10], [11], [11], []]
pathway_exp = [[0, 0.5, 0, 0, 0, 0],
           [-0.41216031, 0.21398791, 0.67943786, 0, 0, 0],
           [0, 0, 0, 0.41216031, 0.21398791, -0.67943786],
           [0, 0, 0, -0.4358056479, 0, 0.6410465853],
           [0, 0, 0, 0.1899306382, 0, 0.4173119146]]
pathway_theo = [[-0.000006496415069, 0.5, 0.000006634790866, 0, 0, 0],
                [-0.42222215, 0.21398549, 0.65800129, 0, 0, 0],
                [0, 0, 0, 0.42221658, 0.21400477, -0.65799561],
                [0, 0, 0, -0.4401418105, 0.000001154439359, 0.6215919761],
                [0, 0, 0, 0.2032369749, 0.00000001405307422, 0.4128421789]]
pathway_new = [[-0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
               [-0.3, -0.5, -0.5, -0.3, -0.5, 0.5],
               [-0.1, -0.5, -0.5, -0.1, -0.5, 0.5],
               [ 0.1, -0.5, -0.5,  0.1, -0.5, 0.5],
               [ 0.3, -0.5, -0.5,  0.3, -0.5, 0.5]]
pathway_new = []
for x in range(10):
    pathway_new += [[-0.5 + 0.25 * i, -0.5 + 0.1 * x, -0.5, -0.5 + 0.25 * i, -0.5 + 0.1 * x, 0.5] for i in range(5)]
names = [['Γ', '1'], ['2', 'Γ'], ['Γ', '3'], ['Γ', '4'], ['Γ', '5']]
names = [['_', '_'] for path in pathway_new]
geo_file = "../../FHI-aims/Yi/Yi_1_5_D/Results/New_Results/n_4/experimental/geometry.in"
pt.plot_3d_solid_with_path_and_names(geo_file, corners, adjacency, pathway_new, names, save=False)
