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
import math

base_path_1 = "../../FHI-aims/Double_Perovskites/New_Structures/Delta_Beta/AgBi/bands/180/10.0/"
base_path_2 = "../../FHI-aims/Double_Perovskites/New_Structures/Delta_Beta/AgBi/bands/180/10.0_hse_test/"
systems = [base_path_1, base_path_2]

base = "/Users/jackmorgenstein/FHI-aims/Double_Perovskites/New_Structures/Delta_Beta/AgBi/test_ref_frame/OOP/"
systems = ["inverted_2/"]
for system in systems:
    continue
    folder = base + system
    pt.mulliken_plot(f'{folder}settings.in', save=False)
    for i in range(1, 6):
        print(bbo.band_info(f'{folder}', f'band100{i}.out', band_gap=False, spin_splitting=True, verbose=False))
    print()

if False:
    def rotate_points(points, axis, angle):
        """
        Rotate a set of 3D points around a given axis vector by a specific angle.

        Parameters:
        - points: Nx3 numpy array of points
        - axis: 3-element iterable representing the rotation axis
        - angle_degrees: angle in degrees

        Returns:
        - rotated_points: Nx3 numpy array of rotated points
        """
        # Convert angle to radians
        angle = np.radians(angle)

        # Normalize the rotation axis
        axis = np.array(axis, dtype=np.float64)
        axis = axis / np.linalg.norm(axis)

        # Rodrigues' rotation formula components
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])

        # Rotation matrix
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

        # Apply rotation
        rotated_points = points @ R.T  # Transpose since we multiply from the right

        return rotated_points


    Bi_octahedra = [[0, 11.787, 0],
                    [-3.003, 11.787, 0],
                    [0, 8.784, 0],
                    [0, 11.787, 3.135],
                    [0, 11.787, -3.135],
                    [3.003, 11.787, 0],
                    [0, 14.79, 0]]
    Bi_octahedra = np.array([np.array(point) - np.array([0, 11.787, 0]) for point in Bi_octahedra])

    Ag_octahedra = [
        [0, 5.894, 0],
        [-2.891, 5.894, 0],
        [0, 3.003, 0],
        [0, 5.893, -3.055],
        [2.891, 5.894, 0],
        [0, 8.784, 0],
        [0, 5.894, 3.055]
    ]
    Ag_octahedra = np.array([np.array(point) - np.array([0, 5.894, 0]) for point in Ag_octahedra])

    # Define rotation axis and angle
    axis = [1, 1, 0]  # Rotate around z-axis

    print("lattice_vector 20 0 0")
    print("lattice_vector 0 20 0")
    print("lattice_vector 0 0 20")

    rotated = rotate_points(rotate_points(Bi_octahedra, axis=[1, 1, 0], angle=30), axis=[0, 0, 1], angle=15)
    for atom in rotated:
        print(f'atom {atom[0] + 10} {atom[1] + 10} {atom[2]} I')

    rotated = rotate_points(rotate_points(Ag_octahedra, axis=[1, -1, 0], angle=30), axis=[0, 0, 1], angle=-15)
    for atom in rotated:
        print(f'atom {atom[0] - 6 + 3.445 - 2.654 + 10} {atom[1] - 0.919 + 0.885 + 10} {atom[2] + 0.04} I')

# base_path = "../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/"
# systems = ["voxkif/", "selwoz/", "selwuf/", "ijayuq/"]

# base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Ref_frame_test/OOP/"
# systems = ["base_bands/", "0_db_out/", "40_db_out/"]

base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Ref_frame_test/beta_avg_test/"
systems = ["base/", "155/", "165/"]

for system in systems:
    for i in range(1, 6):
        print(f'{system} {i}')
        print(bbo.band_info(f'{base_path}{system}', f'band100{i}.out', band_gap=False, spin_splitting=True, verbose=False))
    # pt.mulliken_plot(f'{base_path}{system}settings.in', save=False)
    # bg.random_angle_info(f'{base_path}{system}/geometry.in', method="flat", bonds=[])

raise ValueError("ve")


def plot_data(cur_x, cur_y, x_label):
    plt.clf()
    A = np.vstack([cur_x, np.ones_like(cur_x)]).T
    m, c = np.linalg.lstsq(A, cur_y, rcond=None)[0]
    y_pred = m * cur_x + c
    ss_res = np.sum((cur_y - y_pred) ** 2)
    ss_tot = np.sum((cur_y - np.mean(cur_y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    plt.scatter(cur_x, cur_y, color='r')
    plt.plot(cur_x, y_pred, 'k', label=f"$R^2$ = {r_squared:.3f}")
    plt.xlabel(x_label)
    plt.ylabel(r'$\Delta E^\pm$')
    plt.legend()
    plt.tight_layout()
    plt.show()


def find_max_delta_e(base_path, num_paths):
    spin_splits = []
    for i in range(1, num_paths + 1):
        spin_splits.append(float(bbo.band_info(base_path, f'band100{i}.out', band_gap=False,
                                               spin_splitting=True, verbose=False)))
    return max(spin_splits)


db_old = []
db_in = []
db_out = []
db_new = []
de = []

# VOXKIF
# path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/VOXKIF/geo_plus.in"
# a = ['1', '5', '2']
# b = ['2', '6', '3']
# c = ['3', '7', '4']
# d = ['4', '8', '1']
# de.append(find_max_delta_e("../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/voxkif/", 4))
# db_old.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], method="old"))
# db_in.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], method="in"))
# db_out.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], method="out"))
# db_new.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], method="new"))

# SELWOZ
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/SELWOZ/geo_plus.in"
a = ['1', '4', '8']
b = ['8', '5', '7']
c = ['7', '3', '6']
d = ['6', '2', '1']
de.append(find_max_delta_e("../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/selwoz/", 4))
db_old.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'3': '1,0,0', '5': '1,0,0'}, method="old"))
db_in.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'3': '1,0,0', '5': '1,0,0'}, method="in"))
db_out.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'3': '1,0,0', '5': '1,0,0'}, method="out"))
db_new.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'3': '1,0,0', '5': '1,0,0'}, method="new"))

# SELWUF
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/SELWUF/geo_plus.in"
a = ['1', '5', '2']
b = ['2', '4', '8']
c = ['8', '6', '7']
d = ['7', '3', '1']
de.append(find_max_delta_e("../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/selwuf/", 4))
db_old.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'4': '0,-1,0', '6': '0,-1,0'}, method="old"))
db_in.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'4': '0,-1,0', '6': '0,-1,0'}, method="in"))
db_out.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'4': '0,-1,0', '6': '0,-1,0'}, method="out"))
db_new.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap={'4': '0,-1,0', '6': '0,-1,0'}, method="new"))

# IJAYUQ
path = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/IJAYUQ/geo_plus.in"
a = ['19', '10', '2']
b = ['2', '6', '213']
c = ['213', '11', '1']
d = ['1', '3', '19']
shiftmap = {'1': '0,0,-1', '6': '0,0,-1', '11': '0,0,-1'}
de.append(find_max_delta_e("../../FHI-aims/Double_Perovskites/Real_Structures/calcs_and_outputs/bands/ijayuq/", 4))
db_old.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap, method="old"))
db_in.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap, method="in"))
db_out.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap, method="out"))
db_new.append(bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap, method="new"))

# note IJAYUQ actually has 2 inequivalent cells, I still need to fix this as of 01/09/2026
# also I'm pretty sure my db_out calculations are wrong, need to do both db_out and then avg/max
if True:
    plot_data(np.array(db_old), np.array(de), r'max($\Delta\beta_\text{old}$)')
    # plot_data(np.array(db_in), np.array(de), r'$\Delta\beta_\text{in}$')
    # plot_data(np.array(db_out), np.array(de), r'$\Delta\beta_\text{out}$')
    plot_data(np.array([(elem[0] + elem[1]) / 2 for elem in db_new]), np.array(de), r'$\Delta\beta_\text{avg}$')
    # plot_data(np.array([(2 * elem[0] + elem[1]) / 3 for elem in db_new]), np.array(de), r'$\Delta\beta_\text{new}$')
    # plot_data(np.array([(3 * elem[0] + elem[1]) / 4 for elem in db_new]), np.array(de), r'$\Delta\beta_\text{new}$')
    # plot_data(np.array([(4 * elem[0] + elem[1]) / 5 for elem in db_new]), np.array(de), r'$\Delta\beta_\text{new}$')
    # plot_data(np.array([(5 * elem[0] + elem[1]) / 6 for elem in db_new]), np.array(de), r'$\Delta\beta_\text{new}$')

raise ValueError("ve")

# IDK:
a = [19, 10, 2]
b = [2, 12, 20]
c = [20, 4, 214]
d = [214, 6, 19]
shiftmap = {4: '-1,0,-1', 6: '-1,0,-1', 20: '-1,0,-1'}
# bg.new_robust_delta_beta(path, bonds=[a, b, c, d], shiftmap=shiftmap)

# RECHECK WHAT DELTA BETA OUT IS BEING COMPUTED IN BasicGeo.py

# Random structures:
a = ['4', '9', '1']
b = ['1', '10', '4_2']
c = ['4_2', '7', '1_2']
d = ['1_2', '8', '4']
shiftmap = {'1': '1,1,0', '10': '0,1,0', '4_2': '0,1,0', '7': '0,1,0', '1_2': '0,1,0'}
base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Random/"
# pt.plot_random_structures([f'{base_path}AgBi_static_Cs/'], title="AgBi")
folder = f'{base_path}Pb_static_Cs/'
all_db_old = []
all_db_in = []
all_db_out = []
all_db_new = []
all_de = []
for i in ["05", "10", "15", "20"]:  # iterate through displacement values
    for j in range(1, 11):  # iterate through 10 generated structures
        band_base = folder + i + "/" + str(j).rjust(2, "0") + "/"
        spin_splittings = []
        for k in range(1, 3):  # iterate through band paths
            spin_splittings.append(float(bbo.band_info(band_base, f'band100%s.out' % k, band_gap=False,
                                                       spin_splitting=True, verbose=False)))
        de = max(spin_splittings)
        if de < 0.1:
            all_de.append(de)
            all_db_old.append(
                bg.new_robust_delta_beta(f'{band_base}geometry.in', bonds=[a, b, c, d], shiftmap=shiftmap, method="old"))
            all_db_in.append(
                bg.new_robust_delta_beta(f'{band_base}geometry.in', bonds=[a, b, c, d], shiftmap=shiftmap, method="in"))
            all_db_out.append(
                bg.new_robust_delta_beta(f'{band_base}geometry.in', bonds=[a, b, c, d], shiftmap=shiftmap, method="out"))
            all_db_new.append(
                bg.new_robust_delta_beta(f'{band_base}geometry.in', bonds=[a, b, c, d], shiftmap=shiftmap, method="new"))

plot_data(np.array(all_db_old), np.array(all_de), r'$\Delta\beta_\text{old}$')
plot_data(np.array(all_db_in), np.array(all_de), r'$\Delta\beta_\text{in}$')
plot_data(np.array(all_db_out), np.array(all_de), r'$\Delta\beta_\text{out}$')
plot_data(np.array([(elem[0] + elem[1]) / 2 for elem in all_db_new]), np.array(all_de), r'$\Delta\beta_\text{new}$')
plot_data(np.array([(2 * elem[0] + elem[1]) / 3 for elem in all_db_new]), np.array(all_de),
          r'$\Delta\beta_\text{new, weighted}$')
plot_data(np.array([(3 * elem[0] + elem[1]) / 4 for elem in all_db_new]), np.array(all_de), r'$\Delta\beta_\text{new}$')
plot_data(np.array([(4 * elem[0] + elem[1]) / 5 for elem in all_db_new]), np.array(all_de), r'$\Delta\beta_\text{new}$')
plot_data(np.array([(5 * elem[0] + elem[1]) / 6 for elem in all_db_new]), np.array(all_de), r'$\Delta\beta_\text{new}$')
