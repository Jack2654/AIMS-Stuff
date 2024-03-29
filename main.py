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

set = "20"
base = f'../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/python_dist_positions/%s/' % set
maurer = True
sig_sq = True
delta_d = True
spin_sp = True
debug = False

for i in [10, 11]:
    break
    res = ""
    folder = f'%s%d/' % (base, i)
    geo = folder + "geometry.in"
    control = folder + "control.in"

    shifts_ag1 = [[0, 0, 0] for x in range(6)]
    corners_ag1 = (20, 18, 14, 17, 22, 26)
    center_ag1 = 2

    shifts_ag2 = [[0, 0, 0] for x in range(6)]
    shifts_ag2[2] = [1, 0, 0]
    shifts_ag2[3] = [0, 1, 0]
    corners_ag2 = (13, 15, 19, 16, 21, 25)
    center_ag2 = 1

    shifts_bi1 = [[0, 0, 0] for x in range(6)]
    shifts_bi1[2] = [1, 0, 0]
    corners_bi1 = (14, 16, 20, 15, 23, 27)
    center_bi1 = 3

    shifts_bi2 = [[0, 0, 0] for x in range(6)]
    shifts_bi2[3] = [0, 1, 0]
    corners_bi2 = (19, 17, 13, 18, 24, 28)
    center_bi2 = 4

    if maurer:
        maurer1_ag1, maurer2_ag1 = bg.maurer_displacement(geo, corners_ag1, shifts_ag1, center_ag1, Ag=True)
        maurer1_ag2, maurer2_ag2 = bg.maurer_displacement(geo, corners_ag2, shifts_ag2, center_ag2, Ag=True)
        maurer1_bi1, maurer2_bi1 = bg.maurer_displacement(geo, corners_bi1, shifts_bi1, center_bi1, Ag=False)
        maurer1_bi2, maurer2_bi2 = bg.maurer_displacement(geo, corners_bi2, shifts_bi2, center_bi2, Ag=False)
        if debug:
            print(f'Maurer for %s' % geo)
            print(maurer1_ag1)
            print(maurer1_ag2)
            print(maurer1_bi1)
            print(maurer1_bi2)
            print(maurer2_ag1)
            print(maurer2_ag2)
            print(maurer2_bi1)
            print(maurer2_bi2)
        else:
            maurer1_ag_avg = 0.5 * (maurer1_ag1 + maurer1_ag2)
            maurer1_bi_avg = 0.5 * (maurer1_bi1 + maurer1_bi2)
            maurer1_avg = 0.5 * (maurer1_ag_avg + maurer1_bi_avg)
            maurer2_ag_avg = 0.5 * (maurer2_ag1 + maurer2_ag2)
            maurer2_bi_avg = 0.5 * (maurer2_bi1 + maurer2_bi2)
            maurer2_avg = 0.5 * (maurer2_ag_avg + maurer2_bi_avg)
            res += f'%0.10f %0.10f %0.10f %0.10f %0.10f ' % (maurer1_ag_avg, maurer1_bi_avg, maurer2_ag_avg,
                                                             maurer2_bi_avg, max(maurer1_avg, maurer2_avg))

    if sig_sq:
        sig_sq_ag1 = bg.sigma_squared(geo, corners_ag1, shifts_ag1, center_ag1)
        sig_sq_ag2 = bg.sigma_squared(geo, corners_ag2, shifts_ag2, center_ag2)
        sig_sq_bi1 = bg.sigma_squared(geo, corners_bi1, shifts_bi1, center_bi1)
        sig_sq_bi2 = bg.sigma_squared(geo, corners_bi2, shifts_bi2, center_bi2)
        if debug:
            print(f'Sigma Squared for %s' % geo)
            print(sig_sq_ag1)
            print(sig_sq_ag2)
            print(sig_sq_bi1)
            print(sig_sq_bi2)
        else:
            sig_sq_ag_avg = 0.5 * (sig_sq_ag1 + sig_sq_ag2)
            sig_sq_bi_avg = 0.5 * (sig_sq_bi1 + sig_sq_bi2)
            sig_sq_avg = 0.5 * (sig_sq_ag_avg + sig_sq_bi_avg)
            res += f'%0.10f %0.10f %0.10f ' % (sig_sq_ag_avg, sig_sq_bi_avg, sig_sq_avg)

    if delta_d:
        delta_d_ag1 = bg.delta_d(geo, corners_ag1, shifts_ag1, center_ag1)
        delta_d_ag2 = bg.delta_d(geo, corners_ag2, shifts_ag2, center_ag2)
        delta_d_bi1 = bg.delta_d(geo, corners_bi1, shifts_bi1, center_bi1)
        delta_d_bi2 = bg.delta_d(geo, corners_bi2, shifts_bi2, center_bi2)
        if debug:
            print(f'Delta d for %s' % geo)
            print(delta_d_ag1)
            print(delta_d_ag2)
            print(delta_d_bi1)
            print(delta_d_bi2)
        else:
            delta_d_ag_avg = 0.5 * (delta_d_ag1 + delta_d_ag2)
            delta_d_bi_avg = 0.5 * (delta_d_bi1 + delta_d_bi2)
            delta_d_avg = 0.5 * (delta_d_ag_avg + delta_d_bi_avg)
            res += f'%0.10f %0.10f %0.10f ' % (delta_d_ag_avg, delta_d_bi_avg, delta_d_avg)

    if spin_sp:
        spin_splits = ""
        for j in range(1, 7):
            spin_splits += bbo.band_info(folder, f'band100%d.out' % j, spin_splitting=True, verbose=False) + " "
        res += spin_splits

    if len(res) > 1:
        print(res.strip())

disturbs = ["5/", "10/", "15/", "20/"]

for dist in disturbs:
    break
    # print(dist)
    temp = folder + dist
    for i in range(10):
        current = temp + str(i) + "/"
        for j in range(9, 10):
            # band = "band100" + str(j) + ".out"
            results = ""
            for k in range(8):
                results += str(bg.angle_info(current + "geometry.in", betas[k], shifts[k], directions[k])[2]) + " "
                # print(f'%.10f' % (res[0]))
            print(results)
            # print(f'%0.10f' % bg.delta_d(current + "geometry.in", center_bi2, corners_sig_bi2, debug=False))
            # print(f'%0.10f' % bg.sigma_squared(current + "geometry.in", center_bi2, corners_sig_bi2, debug=False))
            # bbo.band_info(current, band, band_gap=True, spin_splitting=True, verbose=False)

# pt.correlation_plot()

base = "../../FHI-aims/Double_Perovskites/real-systems/experimental-bs-selwoz/"
settings = base + "settings.in"
# pt.mulliken_plot(settings, debug=False, quiet=False, save=False)
for i in range(6):
    band = "band100" + str(i + 1) + ".out"
    # bbo.band_info(base, band, spin_splitting=True)

base = "../../FHI-aims/French_NMR/TA_theoretical/TA_construction/"
one = base + "one_c2.in"
two = base + "two_c2.in"
three = base + "three_c2.in"
temp = base + "s4.in"

# Geometry.combine(one, two, temp, 27, 37)
# Geometry.combine(temp, three, temp, 64, 36)
# bg.recenter(temp, temp, 1, 0, 0, 0)

# folder = "../../FHI-aims/Yi/S_S_Mixing/R_S/exp_bands_w_restart/"
folder = "../../FHI-aims/Yi/Yi_1_5_D/n_4_6/bands/experimental/"
file = folder + "settings_final.in"
for i in range(1, 5):
    band = "band100" + str(i) + ".out"
    # print(bbo.band_info(folder, band, band_gap=True, spin_splitting=False, verbose=False))
# pt.mulliken_plot(file, save=True)

# ots.hybridization()

# ots.plot_dissociation_curves()

ots.TA_NMR_calculations()
