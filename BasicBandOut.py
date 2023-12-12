# This file compiles basic functions for analyzing band*.out output files
# In particular, the following functions exist here:
#
# -> band_info(filepath, band_gap, spin_splitting, effective_mass, verbose, debug):
#           computes info about a given band segment given the prompted flags for what to compute

# imports:
import BasicFunc as bf
import math
import numpy as np
from effmass import inputs, extrema
import BasicGeo as bg


# functions:
def band_info(folder, band, steps=1, band_gap=False, k0=False, spin_splitting=False,
              effective_mass=False, verbose=True, debug=False):
    results = ""
    lines, kpoints, gamma, max_valence, min_conduction, \
        i_conduction_min, j_conduction_min, i_valence_max, j_valence_max = read_band_out(folder + band)

    # sets up arrays of CBM, VBM, and the required number of points around both points set by the "steps" parameter
    points_about_CBM = [0] * (2 * steps + 1)
    points_about_CBM[steps] = min_conduction
    points_about_VBM = [0] * (2 * steps + 1)
    points_about_VBM[steps] = max_valence
    for x in range(steps):
        if i_conduction_min - x <= 0:
            points_about_CBM[steps + (x + 1)] = lines[i_conduction_min + (x + 1)][2 * j_conduction_min + 1]
            points_about_CBM[steps - (x + 1)] = points_about_CBM[steps + (x + 1)]
        elif i_conduction_min + x >= len(lines) - 1:
            points_about_CBM[steps - (x + 1)] = lines[i_conduction_min - (x + 1)][2 * j_conduction_min + 1]
            points_about_CBM[steps + (x + 1)] = points_about_CBM[steps - (x + 1)]
        else:
            points_about_CBM[steps - (x + 1)] = lines[i_conduction_min - (x + 1)][2 * j_conduction_min + 1]
            points_about_CBM[steps + (x + 1)] = lines[i_conduction_min + (x + 1)][2 * j_conduction_min + 1]

        if i_valence_max - x <= 0:
            points_about_VBM[steps + (x + 1)] = lines[i_valence_max + (x + 1)][2 * j_valence_max + 1]
            points_about_VBM[steps - (x + 1)] = points_about_VBM[steps + (x + 1)]
        elif i_valence_max + x >= len(lines) - 1:
            points_about_VBM[steps - (x + 1)] = lines[i_valence_max - (x + 1)][2 * j_valence_max + 1]
            points_about_VBM[steps + (x + 1)] = points_about_VBM[steps - (x + 1)]
        else:
            points_about_VBM[steps - (x + 1)] = lines[i_valence_max - (x + 1)][2 * j_valence_max + 1]
            points_about_VBM[steps + (x + 1)] = lines[i_valence_max + (x + 1)][2 * j_valence_max + 1]

    # pull out the distance in units of reciprocal lattice vectors of the band path
    if gamma == 1:
        kps = kpoints[-1]
    elif gamma > 1:
        kps = kpoints[0]
    else:
        kps = [abs(kpoints[-1][i] - kpoints[0][i]) for i in range(3)]

    # computes the band length in 1/Angstroms by multiplying the fractional distance in a direction with the length of
    # the reciprocal lattice vector in that direction and then finding the overall length of the resultant vector
    recip_vec = bg.reciprocal_vectors(folder + "geometry.in")
    r_len = list(abs(np.linalg.norm(recip_vec[i])) for i in range(3))
    band_length = math.sqrt(sum((r_len[i] * kps[i]) ** 2 for i in range(3)))

    if debug:
        print('Reciprocal lengths: %.5f %.5f %.5f' % (r_len[0], r_len[1], r_len[2]))
        print('Band path direction (reciprocal lattice units): %.6f %.6f %.6f' % (kps[0], kps[1], kps[2]))
        print('Band length (A^-1): %.10f' % band_length)

    # sets up the arrays of x values corresponding to the points around the CBM and VBM
    x_factor = band_length / (len(lines) - 1)
    x_start_CBM = x_factor * (i_conduction_min - steps)
    x_start_VBM = x_factor * (i_valence_max - steps)
    x_vals_CBM = [x_start_CBM + (x_factor * y) for y in range(2 * steps + 1)]
    x_vals_VBM = [x_start_VBM + (x_factor * y) for y in range(2 * steps + 1)]

    if debug:
        print("CBM, surrounding points, and corresponding x values:")
        print(points_about_CBM)
        print(x_vals_CBM)
        print("VBM, surrounding points, and corresponding x values:")
        print(points_about_VBM)
        print(x_vals_VBM)
        print('Band path direction (reciprocal lattice units): %.6f %.6f %.6f' % (kps[0], kps[1], kps[2]))
        print('Band length (A^-1): %.10f' % band_length)

    if band_gap:
        min_coeff = bf.fit_poly(x_vals_CBM, points_about_CBM, 2)
        min_x = bf.minimize_func(min_coeff)

        if min_x > x_vals_CBM[2] or min_x < x_vals_CBM[0]:
            min_x = x_vals_CBM[1]

        min_fit = min_coeff[0] + min_coeff[1] * min_x + min_coeff[2] * min_x * min_x

        max_coeff = bf.fit_poly(x_vals_VBM, points_about_VBM, 2)
        max_x = bf.maximize_func(max_coeff)

        if max_x > x_vals_VBM[2] or max_x < x_vals_VBM[0]:
            max_x = x_vals_VBM[1]

        max_fit = max_coeff[0] + max_coeff[1] * max_x + max_coeff[2] * max_x * max_x

        if debug:
            print('Minimum fit conduction energy: %.10f' % min_fit)
            print('Maximum fit valence energy: %.10f' % max_fit)
        if verbose:
            print('Band gap: %.10f eV' % (min_fit - max_fit))
        else:
            results += ' %.10f' % (min_fit - max_fit)

    if k0:
        if i_conduction_min < len(kpoints) / 2:
            offset = x_vals_CBM[steps]
        else:
            offset = band_length - x_vals_CBM[steps]
        if verbose:
            print('Momentum offset (k0): %.10f A^-1' % offset)
        else:
            results += ' %.10f' % offset

    if spin_splitting:
        points_above_CBM = [0] * (2 * steps + 1)
        points_above_CBM[steps] = lines[i_conduction_min][2 * j_conduction_min + 3]
        for x in range(steps):
            if i_conduction_min - x <= 0:
                points_above_CBM[steps + (x + 1)] = lines[i_conduction_min + (x + 1)][2 * j_conduction_min + 3]
                points_above_CBM[steps - (x + 1)] = points_about_CBM[steps + (x + 1)]
            elif i_conduction_min + x >= len(lines) - 1:
                points_above_CBM[steps - (x + 1)] = lines[i_conduction_min - (x + 1)][2 * j_conduction_min + 3]
                points_above_CBM[steps + (x + 1)] = points_about_CBM[steps - (x + 1)]
            else:
                points_above_CBM[steps - (x + 1)] = lines[i_conduction_min - (x + 1)][2 * j_conduction_min + 3]
                points_above_CBM[steps + (x + 1)] = lines[i_conduction_min + (x + 1)][2 * j_conduction_min + 3]

        min_coeff = bf.fit_poly(x_vals_CBM, points_about_CBM, 2)
        min_x = bf.minimize_func(min_coeff)

        if min_x > x_vals_CBM[2] or min_x < x_vals_CBM[0]:
            min_x = x_vals_CBM[1]

        min_fit = min_coeff[0] + min_coeff[1] * min_x + min_coeff[2] * min_x * min_x

        above_coeff = bf.fit_poly(x_vals_CBM, points_above_CBM, 2)
        above_fit = above_coeff[0] + above_coeff[1] * min_x + above_coeff[2] * min_x * min_x

        if above_fit < min_fit:
            if debug:
                print("Incorrect spin splittings computed as a result of quadratic approximation.")
                print("Defaulting to value derived explicitly by subtracting output energies:")
                print("Spin splitting: %.10f eV" % (points_above_CBM[steps] - points_about_CBM[steps]))
            else:
                results += ' %.10f' % (points_above_CBM[steps] - points_about_CBM[steps])
        else:
            if debug:
                print('Minimum fit conduction energy: %.10f' % min_fit)
                print(points_about_CBM)
                print(points_above_CBM)
                print("min coeff:")
                print(min_coeff)
                print("above coeff:")
                print(above_coeff)
                print('Value of split band immediately above min conduction: %.10f' % above_fit)
            if verbose:
                print('Spin splitting: %.10f eV' % (above_fit - min_fit))
            else:
                results += ' %.10f' % (above_fit - min_fit)

    if effective_mass:
        h_bar = 1.054571817
        m_0 = 9.109383702
        eV_2_J = 1.602176565
        eff_mass_constants = (100 * (h_bar ** 2)) / m_0
        # constant powers of 10:
        # numerator: h_bar^2 = 10^(-68), G^2 = 10^(20)
        # denominator: m_0 = 10^(-31), eV/J conversion = 10^(-19)
        # cancelling these terms leaves a factor of 100 in the numerator

        del_k = list(float(kps[i]) * steps / (len(kpoints) - 1) for i in range(3))
        eff_mass_del_k = [(r_len[i] * del_k[i]) ** 2 for i in range(3)]
        eff_mass_del_E = eV_2_J * (points_about_CBM[0] - 2 * min_conduction + points_about_CBM[-1])

        if eff_mass_del_E == 0:
            # print("NO change in delta E")
            results += ' %.10f' % 10000
        else:
            eff_mass_dir = [eff_mass_constants * x / eff_mass_del_E for x in eff_mass_del_k]
            eff_mass = eff_mass_constants * np.linalg.norm(eff_mass_del_k) / eff_mass_del_E

            if debug:
                print(del_k)
                print("Constants: " + str(eff_mass_constants))
                print("Delta k: " + str(eff_mass_del_k))
                print("Delta E: " + str(eff_mass_del_E))
                print("eV term: " + str(points_about_CBM[0] - 2 * min_conduction + points_about_CBM[-1]))
            if verbose:
                print('Effective Mass Value: %.10f' % eff_mass)
            else:
                results += ' %.10f' % eff_mass

    if not verbose:
        return results.strip()


def read_band_out(file):
    kpoints = []
    gamma = -1
    max_valence = -100
    i_valence_max = -1
    j_valence_max = -1

    min_conduction = 100
    i_conduction_min = -1
    j_conduction_min = -1

    with open(file, "r") as f:
        lines = f.readlines()

    num_same_min = 0
    num_same_max = 0
    for i, k_point in enumerate(lines):
        full = k_point.split()
        kpoints.append([float(full[1]), float(full[2]), float(full[3])])
        if i == 0 or i == len(lines) - 1:
            if "0.0000000" in full[1] and "0.0000000" in full[2] and "0.0000000" in full[3]:
                gamma = float(full[0])
        points = full[4:]
        if not len(points) % 2 == 0:
            print("error in reading output")
            return
        conduction_num = 0
        for j in range(int(len(points) / 2)):
            occupation = points[2 * j]
            energy = points[2 * j + 1]
            if float(occupation) == 1:
                if float(energy) > max_valence:
                    num_same_max = 0
                    max_valence = float(energy)
                    i_valence_max = i
                    j_valence_max = j
                elif float(energy) == max_valence:
                    num_same_max += 1
            elif float(occupation) == 0:
                conduction_num += 1
                # set conduction_num > 4 for m=3
                # set conduction_num > 0 for anything else
                if float(energy) < min_conduction and conduction_num > 0:
                    num_same_min = 0
                    min_conduction = float(energy)
                    i_conduction_min = i
                    j_conduction_min = j
                elif float(energy) == min_conduction:
                    num_same_min += 1
            else:
                print("Partially filled states detected")

    lines_formatted = []
    for line in lines:
        lines_formatted.append([float(x) for x in line.split()[4:]])

    i_conduction_min += int(num_same_min / 2)
    i_valence_max += int(num_same_max / 2)

    return lines_formatted, kpoints, gamma, max_valence, min_conduction, \
        i_conduction_min, j_conduction_min, i_valence_max, j_valence_max


def eff_mass_package(filepath):
    settings = inputs.Settings(extrema_search_depth=0.025, energy_range=0.75, valence_band=False)
    Aims_data = inputs.DataAims(filepath)
    segments = extrema.generate_segments(settings, Aims_data)
    for ele in segments:
        print(str(ele))
        print(ele.finite_difference_effmass())