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
import matplotlib
import numpy as np
import os
import Methods2024 as M24
import scipy

base_path = "../../FHI-aims/Double_Perovskites/old/AgBi-Perovskites/ideal/out/final_OOP/test_ref_frames/"
# pt.mulliken_plot(f'{base_path}base/settings.in', save=False)
# pt.mulliken_plot(f'{base_path}inverted/settings.in', save=False)


base_path = "../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/"
systems = ["voxkif/", "selwoz/", "selwuf/", "ijayuq/"]
for system in systems:
    continue
    print(system)
    # print(bbo.band_info(f'{base_path}{system}', "band1001.out", band_gap=True, spin_splitting=True, verbose=False))
    # pt.mulliken_plot(f'{base_path}{system}settings_zoomed.in', save=False)
    bg.random_angle_info(f'{base_path}{system}/geometry.in', method="flat", bonds=[])

# VOXKIF
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/VOXKIF/geo_plus.in"
a = [1, 5, 2]
b = [2, 6, 3]
c = [3, 7, 4]
d = [4, 8, 1]
bg.new_robust_delta_beta(path, bonds=[a, b, c, d])

# SELWOZ
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/SELWOZ/geo_plus.in"
a = [1, 4, 8]
b = [8, 5, 7]
c = [7, 3, 6]
d = [6, 2, 1]
bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={3: '1,0,0', 5: '1,0,0'})

# SELWUF
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/SELWUF/geo_plus.in"
a = [1, 5, 2]
b = [2, 4, 8]
c = [8, 6, 7]
d = [7, 3, 1]
bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={4: '0,-1,0', 6: '0,-1,0'})

# IJAYUQ
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/IJAYUQ/geo_plus.in"
a = [19, 10, 2]
b = [2, 6, 213]
c = [213, 11, 1]
d = [1, 3, 19]
shiftmap = {1: '0,0,-1', 6: '0,0,-1', 11: '0,0,-1'}
bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap)

a = [19, 10, 2]
b = [2, 12, 20]
c = [20, 4, 214]
d = [214, 6, 19]
shiftmap = {4: '-1,0,-1', 6: '-1,0,-1', 20: '-1,0,-1'}
bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap)

