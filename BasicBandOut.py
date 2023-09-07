# This file compiles basic functions for analyzing band*.out output files
# In particular, the following functions exist here:
#
# -> band_info(filepath, band_gap, spin_splitting, effective_mass, verbose, debug):
#           computes info about a given band segment given the prompted flags for what to compute

# imports:
import BasicFunc as bf
import math
import matplotlib.pyplot as plt
import numpy as np
from effmass import inputs, extrema
import BasicGeo as bg


# functions:
def band_info(folder, band, band_gap=True, spin_splitting=True, effective_mass=True, verbose=True, debug=False):
    # band length calculation wrong bc outputs in reciprocal lattice units, not A^-1
    kpoints = []
    results = ""
    gamma = -1
    max_valence = -100
    i_valence_max = -1
    j_valence_max = -1

    min_conduction = 100
    i_conduction_min = -1
    j_conduction_min = -1

    with open(folder + band, "r") as f:
        lines = f.readlines()

    for i, k_point in enumerate(lines):
        full = k_point.split()
        kpoints.append([full[1], full[2], full[3]])
        # if debug:
        #     print("k-point " + str(full[0]) + ", coords: " + str(full[1]) + " " + str(full[2]) + " " + str(full[3]))
        if i == 0 or i == len(lines) - 1:
            if "0.0000000" in full[1] and "0.0000000" in full[2] and "0.0000000" in full[3]:
                gamma = float(full[0])
        points = full[4:]
        if not len(points) % 2 == 0:
            print("error in reading output")
            return
        for j in range(int(len(points) / 2)):
            occupation = points[2 * j]
            energy = points[2 * j + 1]
            if float(occupation) == 1:
                if float(energy) > max_valence:
                    max_valence = float(energy)
                    i_valence_max = i
                    j_valence_max = j
            elif float(occupation) == 0:
                if float(energy) < min_conduction:
                    min_conduction = float(energy)
                    i_conduction_min = i
                    j_conduction_min = j
            #else:
            #    print("Partially filled states detected")

    if i_conduction_min == 0:
        min_right = float(lines[i_conduction_min + 1].split()[2 * j_conduction_min + 5])
        min_left = min_right
    elif i_conduction_min == len(lines) - 1:
        min_left = float(lines[i_conduction_min - 1].split()[2 * j_conduction_min + 5])
        min_right = min_left
    else:
        min_left = float(lines[i_conduction_min - 1].split()[2 * j_conduction_min + 5])
        min_right = float(lines[i_conduction_min + 1].split()[2 * j_conduction_min + 5])

    if i_valence_max == 0:
        max_right = float(lines[i_valence_max + 1].split()[2 * j_valence_max + 5])
        max_left = max_right
    elif i_valence_max == len(lines) - 1:
        max_left = float(lines[i_valence_max - 1].split()[2 * j_valence_max + 5])
        max_right = max_left
    else:
        max_left = float(lines[i_valence_max - 1].split()[2 * j_valence_max + 5])
        max_right = float(lines[i_valence_max + 1].split()[2 * j_valence_max + 5])

    if debug:
        print("\nmax left: " + str(max_left))
        print('Valence maxima of %.5f found at k_point %d' % (max_valence, i_valence_max + 1))
        print("max right: " + str(max_right))
        print("min left: " + str(min_left))
        print('Conduction minima of %.5f found at k_point %d' % (min_conduction, i_conduction_min + 1))
        print("min right: " + str(min_right))

    if gamma == 1:
        kps = kpoints[len(kpoints) - 1]
    elif gamma == -1:
        kps = [abs(float(kpoints[-1][i]) - float(kpoints[0][i])) for i in range(3)]
    else:
        kps = kpoints[0]
    recip_vec = bg.reciprocal_vectors(folder + "geometry.in")
    r_len = list(abs(np.linalg.norm(recip_vec[i])) for i in range(3))
    band_length = math.sqrt(sum(r_len[i] ** 2 * float(kps[i]) ** 2 for i in range(3)))
    # print("Band length: " + str(band_length))

    if band_gap:
        x_factor = band_length / (len(lines) - 1)
        x_base = x_factor * i_conduction_min
        min_coeff = bf.fit_poly([x_base - x_factor, x_base, x_base + x_factor], [min_left, min_conduction, min_right],
                                2)
        min_x = bf.minimize_func(min_coeff)
        min_fit = min_coeff[0] + min_coeff[1] * min_x + min_coeff[2] * min_x * min_x
        # problem case when min_left = min_conduction = min_right
        x_factor = band_length / (len(lines) - 1)
        x_base = x_factor * i_valence_max
        max_coeff = bf.fit_poly([x_base - x_factor, x_base, x_base + x_factor], [max_left, max_valence, max_right], 2)
        max_x = bf.maximize_func(max_coeff)
        max_fit = max_coeff[0] + max_coeff[1] * max_x + max_coeff[2] * max_x * max_x

        if min_conduction == min_left and min_conduction == min_right:
            min_fit = min_conduction
            min_x = x_factor * i_conduction_min

        if max_valence == max_left and max_valence == max_right:
            max_fit = max_valence
            max_x = x_factor * i_valence_max

        if verbose:
            print('Band gap: %.10f' % (min_fit - max_fit))
        else:
            results += ' %.10f' % (min_fit - max_fit)

    if spin_splitting:
        above_center = float(lines[i_conduction_min].split()[2 * j_conduction_min + 7])
        if i_conduction_min == 0:
            above_right = float(lines[i_conduction_min + 1].split()[2 * j_conduction_min + 7])
            above_left = above_right
        elif i_conduction_min == len(lines) - 1:
            above_left = float(lines[i_conduction_min - 1].split()[2 * j_conduction_min + 7])
            above_right = above_left
        else:
            above_left = float(lines[i_conduction_min - 1].split()[2 * j_conduction_min + 7])
            above_right = float(lines[i_conduction_min + 1].split()[2 * j_conduction_min + 7])

        x_factor = band_length / (len(lines) - 1)
        x_base = x_factor * i_conduction_min
        min_coeff = bf.fit_poly([x_base - x_factor, x_base, x_base + x_factor], [above_left, above_center, above_right], 2)
        spin_sp = min_coeff[0] + min_coeff[1] * min_x + min_coeff[2] * min_x * min_x - min_fit

        if debug:
            print(min_fit)
            print(spin_sp+min_fit)
        if verbose:
            print('Spin splitting: %.10f' % spin_sp)
        else:
            results += ' %.10f' % spin_sp

    if effective_mass:
        if debug:
            print("\nStarting Effective Mass Calculation using Finite-Forward Different Approximations")
            print('Conduction minima of %.5f found at k_point %d' % (min_conduction, i_conduction_min + 1))
            print('k-point coordinates are %.5f %.5f %.5f' % (float(kpoints[i_conduction_min][0]),
                                                              float(kpoints[i_conduction_min][1]),
                                                              float(kpoints[i_conduction_min][2])))
        h_bar = 1.054571817
        m_0 = 9.109383702
        eV_2_J = 1.602176565
        eff_mass_constants = (100 * (h_bar ** 2)) / m_0   # factor of 100 left over after cancelling out powers of 10

        del_k = list(float(kps[i]) / (len(kpoints) - 1) for i in range(3))
        eff_mass_del_k = sum(r_len[i] ** 2 * del_k[i] ** 2 for i in range(3))
        eff_mass_del_E = eV_2_J * (min_left - 2 * min_conduction + min_right)

        eff_mass = eff_mass_constants * eff_mass_del_k / eff_mass_del_E
        if debug:
            print(del_k)
            print("Constants: " + str(eff_mass_constants))
            print("Delta k: " + str(eff_mass_del_k))
            print("Delta E: " + str(eff_mass_del_E))
        if verbose:
            print('Effective Mass Value: %.10f' % eff_mass)
        else:
            results += ' %.10f' % eff_mass

    if not verbose:
        print(results.strip())


def eff_mass_package(filepath):
    settings = inputs.Settings(extrema_search_depth=0.025, energy_range=0.75, valence_band=False)
    Aims_data = inputs.DataAims(filepath)
    segments = extrema.generate_segments(settings, Aims_data)
    for ele in segments:
        print(str(ele))
        print(ele.finite_difference_effmass())

# file2 = "../../FHI-aims/data_aims/Ge_nsp_aims"
# file3 = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_2_4/experimental"
# bbo.eff_mass(file)
# bbo.eff_mass(file2)
# bbo.eff_mass(file3)
