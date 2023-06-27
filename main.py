import numpy as np
import PlottingTools as pt
import Geometry
import KellerPBE as kp
import BasicGeo as bg
import BasicFunc as bf
import VectorToolkit as vt
import math
import time

# for setting up/ running dissociation curve calculations
setting = 6

if setting == 1:
    file_m4_40 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4_40/"
    file_m4_50 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4_50/"
    control_ts = "../../FHI-aims/KellerPBE/control_files/control_pbe_ts.in"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    addtl = "vdw_damping_sr 1.0660114273175167\n vdw_damping_d "
    files = [file_m4_40, file_m4_50]
    for f in files:
        temp = addtl + str(float(f[-3:-1]))
        kp.write_controls_to_dc(f, control_ts, min_defaults, temp)

if setting == 2:
    poss = ["m4_40/", "m4_50/"]
    base = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_"
    times = []
    for p in poss:
        temp = base + p
        temp_t = time.time()
        kp.run_all_dc(temp)
        times.append(time.time() - temp_t)
    for t in times:
        print("Operation took " + str(t) + " seconds")

if setting == 3:
    poss = ["m4_40/", "m4_50/"]
    base = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_"

    for p in poss:
        temp = base + p
        output_file = temp + "binding_energies.txt"
        kp.binding_energies(temp, output_file)

if setting == 4:
    poss = ["m4_40/", "m4_50/", "m4_30/"]
    base = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_"
    file = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight_ts/binding_energies.txt"

    for p in poss:
        temp = base + p
        output_file = temp + "binding_energies.txt"
        print(kp.mean_absolute_error(output_file, file))

if setting == 5:
    poss = ["m4_40/", "m4_50/"]
    base = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_"

    for p in poss:
        temp = base + p
        output_file = temp + "binding_energies.txt"
        res = kp.determine_all_minima(output_file)
        print(p[-5:])
        for r in res:
            print(r)
        print()
        print()

if setting == 6:
    pbe_tight = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight/"
    pbe_tight_ts = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight_ts/"
    min_s = "../../FHI-aims/KellerPBE/dissociation_curves/min_s/"
    min_s_m4 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4/"
    min_s_m4_40 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4_40/"
    # min_s_m4_50 = "../../FHI-aims/KellerPBE/d_sr_optimization/geos_m4_50/"
    series = [pbe_tight, pbe_tight_ts, min_s, min_s_m4, min_s_m4_40]
    kp.plot_dissociation_curve(series, 36)
