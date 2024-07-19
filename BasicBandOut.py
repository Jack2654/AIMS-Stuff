# This file compiles basic functions for analyzing band*.out output files
# In particular, the following functions exist here:
#
# -> band_info(filepath, band_gap, spin_splitting, effective_mass, verbose, debug):
#           computes info about a given band segment given the prompted flags for what to compute

# imports:
import BasicFunc as bf
import math
import numpy as np
# from effmass import inputs, extrema
import BasicGeo as bg
import matplotlib.pyplot as plt


# functions:
def band_info(folder, band, steps=1, band_gap=True, k0=False, spin_splitting=False,
              effective_mass=False, verbose=True, debug=False, display=False):
    results = ""
    lines, kpoints, gamma, max_valence, min_conduction, \
    i_conduction_min, j_conduction_min, i_valence_max, j_valence_max = read_band_out(folder + band)
    print(kpoints)
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

        if display:
            xrange = np.arange(min(x_vals_CBM), max(x_vals_CBM), (max(x_vals_CBM) - min(x_vals_CBM)) / 50)
            print(points_about_CBM)
            print(points_above_CBM)
            plt.plot(xrange, [min_coeff[0] + min_coeff[1] * x + min_coeff[2] * x * x for x in xrange])
            plt.plot(xrange, [above_coeff[0] + above_coeff[1] * x + above_coeff[2] * x * x for x in xrange])
            plt.scatter(x_vals_CBM, points_about_CBM)
            plt.scatter(x_vals_CBM, points_above_CBM)
            plt.axvline(min_x)
            plt.show()

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
            if verbose:
                # print("NO change in delta E")
                print("Effective Mass Value: inf")
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
            elif float(occupation) == 0 or float(occupation) < 1:
                if float(occupation) > 0:
                    print("Partially filled states detected")
                conduction_num += 1
                # set conduction_num > 4 for m=3 to skip I2 interstitial
                # set conduction_num > 0 for anything else
                if float(energy) < min_conduction and conduction_num > 0:
                    num_same_min = 0
                    min_conduction = float(energy)
                    i_conduction_min = i
                    j_conduction_min = j
                elif float(energy) == min_conduction:
                    num_same_min += 1
            else:
                # continue
                # print(k_point)
                print(occupation)
                print("Partially filled states detected")

    lines_formatted = []
    for line in lines:
        lines_formatted.append([float(x) for x in line.split()[4:]])

    i_conduction_min += int(num_same_min / 2)
    i_valence_max += int(num_same_max / 2)

    return lines_formatted, kpoints, gamma, max_valence, min_conduction, \
           i_conduction_min, j_conduction_min, i_valence_max, j_valence_max


def eff_mass_package(filepath):
    # settings = inputs.Settings(extrema_search_depth=0.025, energy_range=0.75, valence_band=False)
    # Aims_data = inputs.DataAims(filepath)
    # segments = extrema.generate_segments(settings, Aims_data)
    # for ele in segments:
    #     print(str(ele))
    #     print(ele.finite_difference_effmass())
    print(0)


def band_gap(folder, band, valence_offset=0, conduction_offset=0, display=False):
    samples, kpoints, cutoff = [], [], 0
    with open(folder + band, "r") as f:
        for line in f.readlines():
            temp = line.split()
            samples.append(temp[4:])
            kpoints.append([float(x) for x in temp[1:4]])
    states = [[float(sample[x]) for sample in samples] for x in range(len(samples[0]))]
    for i_state in range(len(states) // 2):
        occupation = sum(states[2 * i_state]) / len(states[2 * i_state])
        if occupation == 0.0 and cutoff == 0:
            cutoff = i_state
        elif 0.0 < occupation < 1.0:
            print(f'Partially filled states detected at state %s' % (i_state + 1))
            print(states[2 * i_state])
    shift = max(states[2 * cutoff - 1])
    states = [[val - shift for val in state] for state in states]

    recip_vec = bg.reciprocal_vectors(folder + "geometry.in")
    r_len = [abs(np.linalg.norm(recip_vec[i])) for i in range(3)]
    x_vals = [np.linalg.norm([abs(val) * r_len[i] for i, val in enumerate(kp)]) for kp in kpoints]

    cond_state = states[2 * cutoff + 1 + 2 * conduction_offset]
    val_state = states[2 * cutoff - 1 - 2 * valence_offset]
    min_x, fit_min, min_coeff, min_eff = data_fit(x_vals, cond_state, kpoints, r_len, conduction=True)
    max_x, fit_max, max_coeff, max_eff = data_fit(x_vals, val_state, kpoints, r_len, conduction=False)

    if display:
        plt.plot(x_vals, states[2 * cutoff - 1], 'k')
        plt.plot(x_vals, states[2 * cutoff + 1], 'k')
        step_size = (x_vals[1] - x_vals[0]) / 10
        new_x_range = []
        for x in x_vals:
            new_x_range = new_x_range + [x - i * step_size for i in range(10)]
        plt.scatter(new_x_range, [min_coeff[0] + min_coeff[1] * x + min_coeff[2] * x * x for x in new_x_range],
                 color='orange', s=1)
        plt.scatter(new_x_range, [max_coeff[0] + max_coeff[1] * x + max_coeff[2] * x * x for x in new_x_range],
                 color='orange', s=1)
        plt.scatter(min_x, fit_min, color='purple', s=20)
        plt.scatter(max_x, fit_max, color='purple', s=20)
        plt.scatter(x_vals, states[2 * cutoff - 1 - 2 * valence_offset], color='r', s=5)
        plt.scatter(x_vals, states[2 * cutoff + 1 + 2 * conduction_offset], color='b', s=5)
        plt.ylim([-1, 4])
        plt.ylabel("Energy (eV)")
        plt.axhline()
        plt.show()

    return fit_min - fit_max, min_eff, max_eff


def effective_mass(poly_coeff, kpt, num_kpts, r_len, y_vals):
    h_bar = 1.054571817
    m_0 = 9.109383702
    eV_2_J = 1.602176565
    eff_mass_constants = (100 * (h_bar ** 2)) / m_0
    # constant powers of 10:
    # numerator: h_bar^2 = 10^(-68), G^2 = 10^(20)
    # denominator: m_0 = 10^(-31), eV/J conversion = 10^(-19)
    # cancelling these terms leaves a factor of 100 in the numerator

    # steps = (len(y_vals) - 1) / 2
    # del_k = list(kpt[i] * steps / num_kpts for i in range(3))
    # eff_mass_del_k = [(r_len[i] * del_k[i]) ** 2 for i in range(3)]
    # eff_mass_del_E = eV_2_J * (y_vals[0] - 2 * min(y_vals) + y_vals[-1])
    # if eff_mass_del_E == 0:
    #     print(y_vals)
    # eff_mass = eff_mass_constants * np.linalg.norm(eff_mass_del_k) / eff_mass_del_E

    dE_dk = eV_2_J * poly_coeff[2] * 2
    eff_mass = eff_mass_constants / dE_dk

    return eff_mass


def data_fit(x_vals, data, kpoints, r_len, conduction=True):
    extrema = 10000 if conduction else -10000
    startdex = -1
    enddex = -1
    for index, val in enumerate(data):
        if conduction:
            if val < extrema:
                extrema = val
                startdex = index
                enddex = index
            elif val == extrema:
                enddex = index
        else:
            if val > extrema:
                extrema = val
                startdex = index
                enddex = index
            elif val == extrema:
                enddex = index
    step_size = x_vals[1] - x_vals[0]
    x_range = [x_vals[startdex] - 2 * step_size, x_vals[startdex] - step_size] + x_vals[startdex: enddex + 1] + \
              [x_vals[enddex] + step_size, x_vals[enddex] + 2 * step_size]
    if startdex >= 2:
        y_vals = data[startdex - 2: enddex + 1]
    elif startdex == 1 and enddex >= len(data) - 2:
        y_vals = [data[startdex - 1]] + data[startdex - 1: enddex + 1]
    elif startdex == 1:
        y_vals = [data[enddex + 2]] + data[startdex - 1: enddex + 1]
    elif startdex == 0 and enddex == len(data) - 1:
        y_vals = data[startdex: startdex + 2] + data[startdex: enddex + 1]
    elif startdex == 0 and enddex == len(data) - 2:
        y_vals = [data[enddex + 1]] * 2 + data[startdex: enddex + 1]
    else:
        y_vals = [data[enddex + 2], data[enddex + 1]] + data[startdex: enddex + 1]

    if enddex <= len(data) - 3:
        y_vals = y_vals + data[enddex + 1: enddex + 3]
    elif enddex == len(data) - 2 and startdex <= 1:
        y_vals = y_vals + data[enddex + 1] * 2
    elif enddex == len(data) - 2:
        y_vals = y_vals + [data[enddex + 1]] + [data[startdex - 2]]
    elif enddex == len(data) - 1 and startdex == 0:
        y_vals = y_vals + data[startdex: startdex + 2]
    elif enddex == len(data) - 1 and startdex == 1:
        y_vals = y_vals + [data[startdex - 1]] * 2
    else:
        y_vals = y_vals + [data[startdex - 1], data[startdex - 2]]

    poly_coeff = bf.fit_poly(x_range, y_vals, 2)
    extrema_x = bf.minimize_func(poly_coeff) if conduction else bf.maximize_func(poly_coeff)
    direct_extrema = min(data) if conduction else max(data)
    fit_ext = poly_coeff[0] + poly_coeff[1] * extrema_x + poly_coeff[2] * extrema_x * extrema_x
    plt.scatter(extrema_x, fit_ext, color='g', s=20)
    if conduction and (fit_ext > direct_extrema or fit_ext < direct_extrema - 0.01):
        fit_ext = direct_extrema
        extrema_x = x_vals[data.index(direct_extrema)]
    elif not conduction and (fit_ext < direct_extrema or fit_ext > direct_extrema + 0.01):
        fit_ext = direct_extrema
        extrema_x = x_vals[data.index(direct_extrema)]

    gamma = sum(kpoints[0])
    kpt = kpoints[0]
    if gamma == 0:
        kpt = kpoints[-1]

    eff_mass = 0
    if conduction:
        eff_mass = effective_mass(poly_coeff, kpt, len(kpoints) - 1, r_len, y_vals)

    return extrema_x, fit_ext, poly_coeff, eff_mass
