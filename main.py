import numpy as np
import PlottingTools as pt
import Geometry
import KellerPBE as kp
import BasicGeo as bg
import BasicFunc as bf
import BasicAimsOut as bao
import VectorToolkit as vt
import math
import time
import BayesianOptimization as bo
import FlawedPlotting as fp
from ase.io import read, write

setting = 0
# settings:
# -> 0: do nothing
# -> 1: run a set of s66 structures and generate output
# -> 2: plot a given structure from s66
# -> 3: find next node to evaluate S66 for Bayesian optimization
# -> 4: d442 run of structures
# -> 5: binding energies of d442 calculation
# -> 6: plot dos of n_3 experimental atom projected dos
# -> 7: histogram for D442x10
# -> 8: histogram for S66x21
# -> 9: plot a given structure from d442

if setting == 1:
    poss = [["1.01", "20"], ["1.02", "20"]]
    base = "../../FHI-aims/KellerPBE/S66/blyp/geos_"
    control = "../../FHI-aims/KellerPBE/control_files/control_blyp_ts.in"
    # tight_defaults = "../../FHI-aims/Repository/species_defaults/defaults_2020/tight/"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    kp.s66x21_run(poss, base, control, min_defaults, ignore=False, write_control=True)

if setting == 2:
    pbe_tight = "../../FHI-aims/KellerPBE/S66/dissociation_curves/pbe_tight/"
    pbe_tight_ts = "../../FHI-aims/KellerPBE/S66/dissociation_curves/pbe_tight_ts/"
    min_s = "../../FHI-aims/KellerPBE/S66/dissociation_curves/min_s/"
    blyp_min = "../../FHI-aims/KellerPBE/S66/blyp/geos_min_s/"
    series = [pbe_tight, pbe_tight_ts, min_s, blyp_min]
    for i in range(66):
        kp.plot_dissociation_curve_s66(series, i+1)

if setting == 3:
    file = "../../FHI-aims/KellerPBE/info/d_sr_opt.txt"
    bo.next_node(file)

if setting == 4:
    poss = ["blyp_min_s/"]
    base = "../../FHI-aims/KellerPBE/D442/geo_comparison/"
    commands = [""]
    control = "../../FHI-aims/KellerPBE/control_files/control_blyp.in"
    # control_ts = "../../FHI-aims/KellerPBE/control_files/control_pbe_ts.in"
    pbe_tight_defaults = "../../FHI-aims/Repository/species_defaults/defaults_2020/tight/"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    defaults = [min_defaults]

    kp.d442x10_run(poss, base, commands, control, defaults, ignore=False, write_control=False)

if setting == 5:
    folder_1 = "../../FHI-aims/KellerPBE/D442/geo_comparison/pbe_tight_ts/"
    folder_2 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s/"
    folder_3 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s_ts/"
    folder_4 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s_1.02_50/"
    folders = [folder_1, folder_2, folder_3, folder_4]
    for file in folders:
        temp = file + "binding_energies.txt"
        # kp.binding_energies(file, temp)
        res = kp.determine_all_minima_d442(temp, option=2)
        print(temp)
        for r in res:
            print(r)
        print()

if setting == 6:
    file_base = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_3/experimental_atom_proj/atom_proj_dos_I0"
    all_numbers = [3, 4, 5, 6, 7, 8, 9, 10, 56, 57, 58, 59, 60, 61, 62, 63, 109, 110, 111, 112, 113, 114, 115, 116, 163,
               164, 165, 166, 167, 168, 169, 170]
    interstitial_numbers = [10, 63, 116, 170]
    non_interstitial_num = [3, 4, 5, 6, 7, 8, 9, 56, 57, 58, 59, 60, 61, 62, 109, 110, 111, 112, 113, 114, 115, 163,
               164, 165, 166, 167, 168, 169]
    numbers = non_interstitial_num
    files = []
    suffix = ".dat"
    # big contrib:
    # 10
    # 63
    # 116
    # 170
    for n in numbers:
        if n >= 100:
            files.append(file_base + str(n) + suffix)
        elif n >= 10:
            files.append(file_base + "0" + str(n) + suffix)
        else:
            files.append(file_base + "00" + str(n) + suffix)
    pt.dos_plot(files)

if setting == 7:
    base = "../../FHI-aims/KellerPBE/D442/geo_comparison/"
    base2 = "../../FHI-aims/KellerPBE/D442/sr_opt/"
    poss = ["min_s/", "min_s_ts/", "min_s_m4/", "min_s_1.00_70/", "min_s_095/", "min_s_100/"]
    series = []
    for p in poss:
        series.append(base + p)
    comparison = base + "pbe_tight_ts/"
    kp.plot_d442_histogram([base2 + poss[5]], comparison)

if setting == 8:
    pbe_tight_ts = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight_ts/"
    min_s = "../../FHI-aims/KellerPBE/dissociation_curves/min_s/"
    min_s_ts = "../../FHI-aims/KellerPBE/dissociation_curves/min_s_ts/"
    min_s_opt = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_1.00_70/"
    series = [min_s, min_s_ts, min_s_opt]
    kp.plot_s66_histogram(series, pbe_tight_ts)

if setting == 9:
    base = "../../FHI-aims/KellerPBE/D442/geo_comparison/"
    # pbe_tight = "pbe_tight/"
    blyp_min_s = "blyp_min_s/"
    pbe_tight_ts = "pbe_tight_ts/"
    min_s = "min_s/"
    min_s_d_sr_opt = "min_s_1.00_70/"
    ends = [blyp_min_s, pbe_tight_ts, min_s, min_s_d_sr_opt]
    series = []
    for e in ends:
        series.append(base + e)
    for i in range(1, 442):
        kp.plot_dissociation_curve_d442(series, i)

# file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_2_4/experimental/band1001.out"
# bao.analyze_band_path(file)

#file = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/"
#for i in range(10):
#    # print("Structure: " + str(i))
#    temp = file + "20/" + str(i) + "/band100"
#    temp_arr = []
#    for j in range(1, 5):
#        temp_2 = temp + str(j) + ".out"
#        temp_arr.append(Geometry.bandProcess(temp_2, True, 0))
#    print(" ".join(temp_arr))

# 1D Ethan Project
structure = 3
if structure == 2:
    file = "../../FHI-aims/1D_Ethan/2"
    figure_loc = "../../FHI-aims/1D_Ethan/2.png"
    eshift = -0.4913
    ymin = -1
    ymax = 5
    substate = 0
    cd = {"Bi": "m", "Br": "g", "N": "b", "C": "y", "H": "c"}
    labels = ('$\Gamma$', '$1\|\Gamma$', '$2\|\Gamma$', '$3\|4$', '$\Gamma$', '$5$')
    title = "1D Structure 2"
    equal = True
    debug = False

if structure == 3:
    file = "../../FHI-aims/1D_Ethan/3"
    figure_loc = "../../FHI-aims/1D_Ethan/3.png"
    eshift = -0.4233
    ymin = -1
    ymax = 4
    substate = 0
    cd = {"Bi": "m", "I": "g", "N": "b", "C": "y", "H": "c"}
    labels = ('$\Gamma$', '$1\|\Gamma$', '$2\|\Gamma$', '$3\|4$', '$\Gamma$', '$5$')
    title = "1D Structure 3"
    equal = True
    debug = False


# pt.mulliken_plot(file, filename=figure_loc, energyshift=eshift, ymin=ymin, ymax=ymax, substate=substate,
#           color_dict=cd, labels=labels, title=title, eq=equal, debug=debug)

print(Geometry.bandProcess(file + "/band1001.out", True, 0))
