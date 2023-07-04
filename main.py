import numpy as np
import PlottingTools as pt
import Geometry
import KellerPBE as kp
import BasicGeo as bg
import BasicFunc as bf
import VectorToolkit as vt
import math
import time
import BayesianOptimization as bo

setting = 7
# settings:
# -> 0: do nothing
# -> 1: run a set of structures and generate output
# -> 2: plot a given structure from s66
# -> 3: find next node to evaluate S66 for Bayesian optimization
# -> 4: d442 run of structures
# -> 5: total energies of d442 calculation
# -> 6: temp
# -> 7: plot dos of n_3 experimental atom projected dos

if setting == 1:
    poss = [[0.91870443, 45]]
    base = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_"
    control_ts = "../../FHI-aims/KellerPBE/control_files/control_pbe_ts.in"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    kp.s66x21_run(poss, base, control_ts, min_defaults, ignore=False)

if setting == 2:
    pbe_tight = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight/"
    pbe_tight_ts = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight_ts/"
    min_s = "../../FHI-aims/KellerPBE/dissociation_curves/min_s/"
    min_s_m4 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_1.03_50/"
    min_s_m4_40 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_1.00_70/"
    # min_s_m4_50 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4_50/"
    series = [pbe_tight, pbe_tight_ts, min_s, min_s_m4, min_s_m4_40]
    kp.plot_dissociation_curve(series, 59)

if setting == 3:
    bo.next_node("../../FHI-aims/KellerPBE/info/d_sr_opt.txt")

if setting == 4:
    folder = "../../FHI-aims/KellerPBE/D442/geo_comparison/pbe_tight_ts/"
    control = "../../FHI-aims/KellerPBE/control_files/control_pbe_ts.in"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    pbe_tight_defaults = "../../FHI-aims/Repository/species_defaults/defaults_2020/tight/"
    kp.write_controls_to_dc(folder, control, pbe_tight_defaults)
    start_time = time.time()
    kp.run_all_dc(folder, ignore=False)
    print("Whole calc took " + str(time.time() - start_time) + " seconds")

if setting == 5:
    folder_1 = "../../FHI-aims/KellerPBE/D442/geo_comparison/pbe_tight_ts/"
    folder_2 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s/"
    folder_3 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s_ts/"
    folder_4 = "../../FHI-aims/KellerPBE/D442/geo_comparison/min_s_1.02_50/"
    folders = [folder_1, folder_2, folder_3, folder_4]
    for file in folders:
        temp = file + "total_energies.txt"
        kp.all_energies(file, temp)

if setting == 6:
    poss = ["min_s/", "min_s_1.02_50/", "min_s_ts/", "pbe_tight_ts/"]
    base = "../../FHI-aims/KellerPBE/D442/geo_comparison/"
    count = 0
    for file in poss:
        cur_folder = base + file
        be_file = cur_folder + "total_energies.txt"
        res = kp.determine_all_minima_d442(be_file)
        print()
        print(file)
        for r in res:
            if r < 0.5 or r > 1.5:
                print(r)
                count += 1
    print(count)

if setting == 7:
    file_base = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_3/experimental_atom_proj/atom_proj_dos_I0"
    numbers = [3, 4, 5, 6, 7, 8, 9, 10, 56, 57, 58, 59, 60, 61, 62, 63, 109, 110, 111, 112, 113, 114, 115, 116, 163, 164, 165, 166, 167, 168, 169, 170]
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











