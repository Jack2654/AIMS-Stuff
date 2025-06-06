# file for 1 time use scripts

# imports
import BasicGeo as bg
import os
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
import math
import BasicAimsOut as bao
from scipy import interpolate
import numpy as np
from scipy.optimize import curve_fit


# code


# standardizes and moves 2D ideal disturbed models:
def move_and_standard(readfile, writefile, cutoff):
    lv = bg.lattice_vectors(readfile)
    at = bg.atoms(readfile)
    ret_at = []
    for atom in at:
        temp = atom.split()
        temp_loc = [float(x) for x in temp[1:4]]
        for i in range(3):
            if temp_loc[i] > cutoff:
                temp_loc[0] -= lv[i][0]
                temp_loc[1] -= lv[i][1]
                temp_loc[2] -= lv[i][2]

        temp_ret = " ".join([str(temp_loc[0]), str(temp_loc[1]), str(temp_loc[2])])
        temp_ret = temp[0] + " " + temp_ret + " " + temp[-1] + "\n"
        ret_at.append(temp_ret)

    lattice_vectors = []
    for lat in lv:
        temp_lat = [str(x) for x in lat]
        temp = "lattice_vector " + " ".join(temp_lat) + "\n"
        lattice_vectors.append(temp)
    with open(writefile, "w") as f:
        f.writelines(lattice_vectors)
        f.writelines(ret_at)


def add_species_default_command(folder):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    for file in all_dir:
        if "read" in file or ".swp" in file:
            continue
        temp = folder + file
        lines = []
        with open(temp, "r") as f:
            for ln in f:
                lines.append(ln)
        with open(temp, "w") as f:
            for line in lines:
                if "species_default_type" in line:
                    f.write("  species_default_type minimal+s\n")
                elif "  e.g." in line:
                    f.write(line)
                    f.write('#  The "species_default_type minimal+s" keyword in this file enables energy, force,\n'
                            + '#  and stress corrections for this given species. See manual for more info.\n')
                else:
                    f.write(line)


def me_331_lab_1_plots():
    pressures = [413.8857256, 2778.947015, 2778.947015, 413.8857256]
    volumes = [0.0002111990915, 0.0002037329002, 0.0002410638566, 0.0002468708943,
               0.0002070512074, 0.0001995850161, 0.000239404703, 0.0002443821639,
               0.0002070512074, 0.0002004145929, 0.0002369159726, 0.0002427230103]
    x_offset = [0.000002, 0.000002, -0.000008, -0.000008]
    y_offset = [100, -130, -130, 100]
    labels = ["Point A", "Point B", "Point C", "Point D"]

    # Cycle 1:
    plt.scatter(volumes[0:4], pressures)
    for i in range(4):
        plt.text(volumes[i] + x_offset[i], pressures[i] + y_offset[i], labels[i])
    plt.plot(volumes[0:2], pressures[0:2], 'b')
    plt.plot(volumes[1:3], pressures[1:3], 'b')
    plt.plot(volumes[2:4], pressures[2:4], 'b')
    plt.plot([volumes[3]] + [volumes[0]], [pressures[3]] + [pressures[0]], 'b')
    plt.xlabel("Volume (m^3)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Part B, Group E, Cycle 1")
    plt.show()

    # Cycle 2:
    plt.scatter(volumes[4:8], pressures)
    for i in range(4, 8):
        plt.text(volumes[i] + x_offset[i - 4], pressures[i - 4] + y_offset[i - 4], labels[i - 4])
    plt.plot(volumes[4:6], pressures[0:2], 'b')
    plt.plot(volumes[5:7], pressures[1:3], 'b')
    plt.plot(volumes[6:8], pressures[2:4], 'b')
    plt.plot([volumes[7]] + [volumes[4]], [pressures[3]] + [pressures[0]], 'b')
    plt.xlabel("Volume (m^3)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Part B, Group E, Cycle 2")
    plt.show()

    # Cycle 3:
    plt.scatter(volumes[8:12], pressures)
    for i in range(8, 12):
        plt.text(volumes[i] + x_offset[i - 8], pressures[i - 8] + y_offset[i - 8], labels[i - 8])
    plt.plot(volumes[8:10], pressures[0:2], 'b')
    plt.plot(volumes[9:11], pressures[1:3], 'b')
    plt.plot(volumes[10:12], pressures[2:4], 'b')
    plt.plot([volumes[11]] + [volumes[8]], [pressures[3]] + [pressures[0]], 'b')
    plt.xlabel("Volume (m^3)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Part B, Group E, Cycle 3")
    plt.show()


def me344_lab5_plots(readfile):
    with open(readfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        time = []
        var_1 = []
        var_2 = []
        var_3 = []
        var_4 = []
        var_5 = []
        var_6 = []
        var_7 = []
        for row in csv_reader:
            time.append(float(row[0]))
            var_1.append(float(row[1]))
            var_2.append(float(row[2]))
            var_3.append(float(row[3]))
            var_4.append(float(row[4]))
            var_5.append(float(row[5]))
            var_6.append(float(row[6]))
            var_7.append(float(row[7]))

    plt.ylim([-10, 250])
    plt.xlim([0, 10])
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (mm)")
    plt.title("Spring C")
    plt.plot(time, var_2)
    plt.plot(time, var_4)
    plt.plot(time, var_3)
    plt.legend(["Road Position", "Middle Position", "Top Position"])
    plt.show()

    print(max(var_3))
    print(min(var_3))
    print(max(var_3) - min(var_3))
    print(max(var_2) - min(var_2))

    if False:
        plt.plot(time, var_5)
        plt.ylim([0, 210])
        plt.xlim([0, 10])
        plt.title("var 5")
        plt.show()

        plt.plot(time, var_6)
        plt.ylim([0, 210])
        plt.xlim([0, 10])
        plt.title("var 6")
        plt.show()

        plt.plot(time, var_7)
        plt.ylim([0, 210])
        plt.xlim([0, 10])
        plt.title("var 7")
        plt.show()


def me321_team_proj_2(file, color):
    with open(file) as f:
        lines = f.readlines()
    reformatted_lines = []
    for ln in lines:
        temp = ln.split()
        temp = [x.replace('\x00', "") for x in temp]
        if len(temp) > 1:
            reformatted_lines.append(temp)
    actual_data = []
    for ln in reformatted_lines:
        if len(ln) == 4:
            if len(ln[0]) > 1 and len(ln[1]) > 1 and len(ln[2]) > 1 and len(ln[3]) > 1:
                actual_data.append([float(x) for x in ln])
    time = [elem[0] for elem in actual_data]
    force = [elem[1] for elem in actual_data]
    displ = [elem[2] for elem in actual_data]
    max = 0
    for f in force:
        if f > max:
            max = f
        if f < max and f > 10:
            break
    print(max)
    start_index = 0
    for i in range(len(force)):
        if force[i] > 1:
            start_index = i - 1
            break
    time = time[start_index:]
    displ = displ[start_index:]
    time = [x - time[0] for x in time]
    displ = [x - displ[0] for x in displ]
    force = force[start_index:]
    plt.plot(displ, force, color, linestyle="-")
    # plt.axhline(max, color=color)


def make_maurer(base):
    files = ["00/", "05/", "10/", "15/", "20/", "25/", "30/", "35/", "40/", "45/", "50/"]
    with open(base + "geometry.in", "r") as f:
        lv = []
        at = []
        for ln in f.readlines():
            if "atom" in ln:
                at.append(ln)
            if "lattice" in ln:
                lv.append(ln)
    for i in range(11):
        temp = i * 0.05
        disp = math.sqrt(0.5 * temp * temp)
        with open(base + files[i] + "geometry.in", "w") as f:
            f.writelines(lv)
            for ln in at:
                if "I" in ln:
                    f.write(ln)
                elif "Cs" in ln:
                    f.write(ln)
                else:
                    temp_at = [float(x) for x in ln.split()[1:4]]
                    temp_at[0] += disp
                    temp_at[1] += disp
                    f.write(f'atom %0.8f %0.8f 0.0 %s\n' % (temp_at[0], temp_at[1], ln.split()[-1]))


def allocate_plotting_files(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    print(all_dir)


def TMS_work():
    folder = "../../FHI-aims/French_NMR/TMS_reference/"
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    for direct in all_dir:
        if os.path.isdir(folder + direct) and "J" not in direct:
            new_direct = " ../../FHI-aims/French_NMR/TMS_full/" + direct
            command = "mkdir " + new_direct
            command2 = "cp " + folder + "u-cc-pVTZ/geometry.in " + new_direct + "/geometry.in"
            command3 = "cp " + folder + direct + "/control.in " + new_direct + "/control.in"
            # os.system(command)
            os.system(command2)
            os.system(command3)


def TMS_work_2():
    folder = "../../FHI-aims/French_NMR/TMS_non_relativistic/"
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    for direct in all_dir:
        if os.path.isdir(folder + direct):
            with open(folder + direct + "/control.in", "r") as f:
                lns = f.readlines()
            lns_to_be_written = []
            lns_to_be_written.append("# Tetramethylsilane - Si(CH3)4\n")
            lns_to_be_written.append("check_stacksize         .false.\n")
            lns_to_be_written.append("xc                      pbe\n")
            lns_to_be_written.append("vdw_correction_hirshfeld\n")
            lns_to_be_written.append("spin                    none\n")
            lns_to_be_written.append("relativistic            none\n")
            lns_to_be_written.append("basis_threshold         0.0\n")
            lns_to_be_written.append("magnetic_response       s\n\n")
            lns_to_be_written.append("# " + direct + "\n")
            species_found = False
            for ln in lns:
                if "species" in ln and not species_found:
                    species_found = True
                if species_found:
                    lns_to_be_written.append(ln)
            with open(folder + direct + "/control.in", "w") as f:
                f.writelines(lns_to_be_written)


def TMS_NMR_shielding_plot(data):
    with open(data, "r") as f:
        basis_type = []
        non_rel_norm = []
        rel_norm = []
        set_size = []
        for line in f:
            temp = line.split()
            basis_type.append(temp[0])
            non_rel_norm.append([float(x) for x in temp[5:8]])
            rel_norm.append([float(x) for x in temp[11:14]])
            set_size.append(int(temp[14]))

    index_dict = {"FHI-aims09": 0, "NAO-VCC-nZ": 1, "pcS-n": 2, "aug-pcS-n": 3,
                  "cc-pVnZ": 4, "aug-cc-pVnZ": 5, "cc-pCVnZ": 6, "aug-cc-pCVnZ": 7}

    color_dict = {"FHI-aims09": 'k', "NAO-VCC-nZ": 'b', "pcS-n": 'saddlebrown', "aug-pcS-n": 'c',
                  "cc-pVnZ": 'y', "aug-cc-pVnZ": 'g', "cc-pCVnZ": 'r', "aug-cc-pCVnZ": 'pink'}

    full_data = [[[], [], [], [], [], [], [], []] for x in range(8)]
    for i in range(len(basis_type)):
        if basis_type[i] in index_dict.keys():
            full_data[index_dict.get(basis_type[i])][0] = basis_type[i]
            full_data[index_dict.get(basis_type[i])][1].append(set_size[i])
            full_data[index_dict.get(basis_type[i])][2].append(non_rel_norm[i][0])  # Si non_rel
            full_data[index_dict.get(basis_type[i])][3].append(non_rel_norm[i][1])  # C non_rel
            full_data[index_dict.get(basis_type[i])][4].append(non_rel_norm[i][2])  # H non_rel
            full_data[index_dict.get(basis_type[i])][5].append(rel_norm[i][0])  # Si rel
            full_data[index_dict.get(basis_type[i])][6].append(rel_norm[i][1])  # C rel
            full_data[index_dict.get(basis_type[i])][7].append(rel_norm[i][2])  # H rel

    titles = ["Si Non-relativistic", "C Non-relativistic", "H Non-relativistic",
              "Si Relativistic", "C Relativistic", "H Relativistic"]
    bounds = [[0.001, 100], [0.006, 20], [0.0001, 3],
              [20, 300], [2, 30], [1.2, 4]]
    for i, cur_title in enumerate(titles):
        for subset in full_data:
            plt.plot(subset[1], subset[i + 2], color=color_dict.get(subset[0]), label=subset[0])
            plt.scatter(subset[1], subset[i + 2], color=color_dict.get(subset[0]), label="_none")
        plt.legend()
        plt.ylim(bounds[i])
        plt.xlabel("Basis Set Size")
        plt.ylabel("error (ppm)")
        plt.semilogy()
        plt.title(cur_title)
        name = "/".join(data.split("/")[:-1]) + "/" + cur_title.replace(" ", "_") + ".png"
        plt.savefig(name, dpi=1000, bbox_inches='tight', format="png")
        plt.clf()


def TA_crystal(additive):
    geo_file = "../../FHI-aims/French_NMR/TA_crystal/geo_full_crystal.in"
    write_file = "../../FHI-aims/French_NMR/TA_crystal/new_geo" + str(additive) + ".in"
    N = []
    H = []
    C = []
    lv = []
    with open(geo_file, "r") as f:
        for ln in f:
            if "H" in ln:
                H.append(ln)
            elif "C" in ln:
                C.append(ln)
            elif "N" in ln:
                N.append(ln)
            elif "lattice_vector" in ln:
                lv.append(ln)
            else:
                print(ln.strip())

    with open(write_file, "w") as f:
        f.writelines(lv)
        for i, ln in enumerate(N):
            if (i + additive) % 4 == 0:
                f.write(ln)
        for i, ln in enumerate(C):
            if (i + additive) % 4 == 0:
                f.write(ln)
        for i, ln in enumerate(H):
            if (i + additive) % 4 == 0:
                f.write(ln)


# computes hybridization of m=4, n=6 system
def hybridization():
    folder = "../../FHI-aims/Yi/Yi_1_5_D/n_4_6/bands/experimental/"
    k_point = []
    found = False
    with open(folder + " bandmlk1001.out", "r") as f:
        for ln in f:
            if "k point number:   101" in ln:
                found = True
            if found and ("4952" in ln or "4953" in ln):
                k_point.append(ln)
    with open(folder + "hybrid_check.out", "w") as f:
        f.writelines(k_point)


def set_up_dissociation_curves(folder, xc, distances):
    for distance in distances:
        current = folder + str(distance)
        os.mkdir(current)
        default_path = "~/FHI-aims/Repositories/FHIaims/species_defaults/defaults_2020/tight/05_B_default"
        with open(current + "/geometry.in", "w") as f:
            f.write("atom 0 0 0 B\n")
            f.write("initial_moment 1\n")
            f.write("atom 0 0 " + str(distance) + " B\n")
            f.write("initial_moment -1")
        with open(current + "/control.in", "w") as f:
            if xc == "pbe_ts":
                print("PBE+TS detected")
                f.write("xc pbe\n")
                f.write("vdw_correction_hirshfeld_alkali\n")
            else:
                f.write("xc " + xc + "\n")
            f.write("relativistic atomic_zora scalar\n")
            f.write("spin collinear\n")
            f.write("check_stacksize .false.\n")
            f.write("sc_accuracy_rho 0.0001\n")
        os.system("cat " + default_path + " >> " + current + "/control.in")


def plot_dissociation_curves():
    distances = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.85, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0]
    folder = "../../Documents/23-24/ME511/In_Class/Ex_1/"
    pw_lda = folder + "pw_lda/"
    pbe = folder + "pbe/"
    pbe_ts = folder + "pbe_ts/"

    # reference energies
    ref_pw_lda = bao.find_total_energy(pw_lda + "free_x_fill/aims.out")
    ref_pbe = bao.find_total_energy(pbe + "free_x_fill/aims.out")
    ref_pbe_ts = bao.find_total_energy(pbe_ts + "free_x_fill/aims.out")
    print(ref_pw_lda, ref_pbe, ref_pbe_ts)

    # total energies
    pw_lda_energies = []
    pbe_energies = []
    pbe_ts_energies = []
    for distance in distances:
        pw_lda_energies.append(float(bao.find_total_energy(pw_lda + "curve/" + str(distance) + "/aims.out")))
        pbe_energies.append(float(bao.find_total_energy(pbe + "curve/" + str(distance) + "/aims.out")))
        pbe_ts_energies.append(float(bao.find_total_energy(pbe_ts + "curve/" + str(distance) + "/aims.out")))

    # binding energies
    binding_pw_lda = [x - 2 * float(ref_pw_lda) for x in pw_lda_energies]
    binding_pbe = [x - 2 * float(ref_pbe) for x in pbe_energies]
    binding_pbe_ts = [x - 2 * float(ref_pbe_ts) for x in pbe_ts_energies]
    print("PBE:")
    print(binding_pbe)
    print("PBE+TS:")
    print(binding_pbe_ts)

    # plotting 3 xc
    spline_pw_lda = interpolate.CubicSpline(distances, binding_pw_lda)
    spline_pbe = interpolate.CubicSpline(distances, binding_pbe)
    spline_pbe_ts = interpolate.CubicSpline(distances, binding_pbe_ts)
    xs = np.arange(1.0, 8.0, 0.01)
    plt.scatter(distances, binding_pw_lda, label='_nolegend_', color='r')
    plt.plot(xs, spline_pw_lda(xs), color='r')
    plt.scatter(distances, binding_pbe, label='_nolegend_', color='b')
    plt.plot(xs, spline_pbe(xs), color='b')
    plt.scatter(distances, binding_pbe, label='_nolegend_', color='k')
    plt.plot(xs, spline_pbe_ts(xs), color='k')
    plt.xlabel("Intermolecular Distance (Angstroms)")
    plt.ylabel("Binding Energy (eV)")
    plt.legend(["pw-lda", "pbe", "pbe+ts"])
    plt.show()

    # plotting pw-lda / LJ / Morse
    epsilon = 3.16612
    sigma = 1.6 / (2 ** (1 / 6))
    alpha = 1.55
    LJ_pw_lda = [4 * epsilon * (((sigma / x)) ** 12 - ((sigma / x) ** 6)) for x in xs]
    morse_pw_lda = [epsilon * ((1 - (math.e ** (alpha * (1.6 - x)))) ** 2 - 1) for x in xs]
    plt.plot(xs, spline_pw_lda(xs))
    plt.plot(xs, LJ_pw_lda)
    plt.plot(xs, morse_pw_lda)
    plt.xlabel("Intermolecular Distance (Angstroms)")
    plt.ylabel("Binding Energy (eV)")
    plt.legend(["pw-lda", "LJ", "Morse"])
    plt.ylim([-4, 8])
    # plt.show()
    plt.clf()

    # plotting pbe / LJ / Morse
    epsilon = 2.59353
    sigma = 1.6 / (2 ** (1 / 6))
    alpha = 1.7
    LJ_pw_lda = [4 * epsilon * (((sigma / x)) ** 12 - ((sigma / x) ** 6)) for x in xs]
    morse_pw_lda = [epsilon * ((1 - (math.e ** (alpha * (1.6 - x)))) ** 2 - 1) for x in xs]
    plt.plot(xs, spline_pbe(xs))
    plt.plot(xs, LJ_pw_lda)
    plt.plot(xs, morse_pw_lda)
    plt.xlabel("Intermolecular Distance (Angstroms)")
    plt.ylabel("Binding Energy (eV)")
    plt.legend(["pbe (+ts)", "LJ", "Morse"])
    plt.ylim([-3, 8])
    plt.show()


# sets up input files for an NMR calculation
def TA_NMR_calculations():
    base = "../../FHI-aims/French_NMR/step3/NMR/"
    TA_folder = base + "I3/"
    TA_geometry = TA_folder + "geometry.in"
    TA_control = TA_folder + "control.in"
    TA_submit = TA_folder + "submit.sh"
    for i in range(63):
        cur_dir = TA_folder + str(i) + "/"
        if i < 10:
            cur_dir = TA_folder + "0" + str(i) + "/"
        cur_geo = cur_dir + "geometry.in"
        cur_sub = cur_dir + "submit.sh"
        command = f'cp %s %scontrol.in' % (TA_control, cur_dir)
        if not os.path.exists(cur_dir):
            os.system("mkdir " + cur_dir)
        os.system(command)

        with open(TA_geometry, "r") as f:
            ln = f.readlines()
        with open(cur_geo, "w") as f:
            count = 0
            for line in ln:
                if "H" in line and count == i:
                    temp = line.split()[:-1]
                    temp = "\t".join(temp)
                    temp += "\tH_NAO"
                    f.write(temp.strip() + "\n")
                    f.write("magnetic_response\n")
                else:
                    f.write(line.strip() + "\n")
                if "H" in line:
                    count += 1

        with open(TA_submit, "r") as f:
            ln = f.readlines()
        with open(cur_sub, "w") as f:
            for line in ln:
                if "name" in line:
                    f.write("#SBATCH --job-name=I3_N_" + str(i) + "\n")
                else:
                    f.write(line)


# removes MPI errors after "Have a nice day" message in aims.out file
def aims():
    folder = "../../FHI-aims/Yi/S_S_Mixing/R_S/exp_bands_w_restart/"
    file = folder + "aims_unaltered.out"
    with open(file, "r") as f:
        lines = f.readlines()
    file = folder + "aims.out"
    with open(file, "w") as f:
        done = False
        for line in lines:
            if done:
                f.write(line)
                break
            if "Have a nice day" in line:
                done = True
            f.write(line)


def plots_for_ME511_ex_2():
    folder = "../../Documents/23-24/ME511/In_Class/Ex_2/calcs/total_energies/"
    constants = []
    energies = []
    for count, i in enumerate(["bcc", "fcc", "diamond"]):
        constants.append([])
        energies.append([])
        all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder + i + "/"))]
        all_dir.sort()
        for dir in all_dir:
            if os.path.isdir(folder + i + "/" + dir):
                constants[count].append(float(dir))
                if i == "diamond":
                    energies[count].append(float(bao.find_total_energy(folder + i + "/" + dir + "/aims.out")) / 2)
                else:
                    energies[count].append(float(bao.find_total_energy(folder + i + "/" + dir + "/aims.out")))
    print(energies)
    constants = [[y ** 3 for y in x] for x in constants]
    constants[2] = [x / 2 for x in constants[2]]
    energies = [[y + 7868.877654655 for y in x] for x in energies]
    print(constants)
    for x in range(3):
        plt.plot(constants[x], energies[x])
    plt.legend(["bcc", "fcc", "diamond"])
    plt.xlabel("Volume per Atom (Angstroms^3)")
    plt.ylabel("Cohesive energy per atom (eV)")
    plt.show()


def more_plots_ME511_ex_2():
    folder = "../../Documents/23-24/ME511/In_Class/Ex_2/calcs/murn/"
    for y in [f'%sfit_%s.dat' % (folder, x) for x in ["bcc", "fcc", "diamond"]]:
        volumes = []
        energies = []
        with open(y, "r") as f:
            for line in f:
                if "#" not in line:
                    temp = line.split()
                    volumes.append(float(temp[0]))
                    energies.append(float(temp[1]) + 7868.877654655)
        plt.plot(volumes, energies)
        plt.legend(["bcc", "fcc", "diamond"])
        plt.xlabel("Volume per Atom (Angstroms^3)")
        plt.ylabel("Cohesive energy per atom (eV)")
    plt.show()


def fix_setting_colors(folder):
    colors = []
    with open(folder + "settings_final.in", "r") as f:
        for line in f.readlines():
            if "color_dict" in line:
                colors = "color_dict      Pb1 pink Pb2 purple I1 lightgreen I2 g I olive " + " ".join(
                    line.split()[5:]) + "\n"
    with open(folder + "settings_plane.in", "r") as f:
        lines = f.readlines()
    with open(folder + "settings_new_plane.in", "w") as f:
        for line in lines:
            if "color" not in line:
                f.write(line)
            else:
                f.write(colors)


def ethanol_j_plots():
    # single ethanol per calculation
    # eth_1_shieldings = [30.29, 30.03, 30.04, 27.19, 27.19, 30.72]
    # eth_2_shieldings = [30.53, 29.90, 30.04, 27.28, 27.75, 31.26]

    # two ethanol in one geometry.in file
    # eth_1_shieldings = [30.19812, 29.91503, 29.5112, 27.04088, 27.12101, 30.23023]
    eth_2_shieldings = [30.19076, 29.93676, 30.0316, 27.23627, 27.24359, 26.76537]

    dmf_shieldings = [28.52817, 27.8061, 28.52855, 29.1626, 29.16272, 26.45794, 22.81182]
    eth_1_shieldings = dmf_shieldings

    eth_1_shieldings = [x - 31.1846 for x in eth_1_shieldings]
    eth_2_shieldings = [x - 31.1846 for x in eth_2_shieldings]
    plt.hist(eth_1_shieldings, bins=[-9 + 0.15 * x for x in range(100)], density=True, rwidth=0.5)
    plt.xlim([-9, 0])
    plt.xlabel("ppm")
    plt.show()
    plt.clf()
    plt.hist(eth_2_shieldings, bins=[-9 + 0.15 * x for x in range(100)], density=True, rwidth=0.5)
    plt.xlim([-5, 0])
    plt.xlabel("ppm")
    plt.show()
    plt.clf()

    eth_1_j = [21.4764, 21.3296, 21.3291, 21.5257]
    eth_2_j = [21.4216, 21.3043, 21.3439, 22.0600]
    plt.hist(eth_1_j, bins=[20 + 0.15 * x for x in range(100)], density=True, rwidth=0.5)
    plt.xlim([21, 22.5])
    plt.xlabel("Hz")
    plt.show()
    plt.clf()
    plt.hist(eth_2_j, bins=[20 + 0.15 * x for x in range(100)], density=True, rwidth=0.5)
    plt.xlim([21, 22.5])
    plt.xlabel("Hz")
    plt.show()
    plt.clf()


def add_magnetic_response(file, write):
    with open(file, "r") as f:
        lines = f.readlines()
    with open(write, "w") as f:
        for line in lines:
            f.write(line)
            if "atom" in line and "H" in line:
                f.write("magnetic_response\n")


# creates relaxations for many FAI possibilities
def make_FAI_relaxations():
    folder = "../../FHI-aims/French_NMR/FA/FAI/relax_many/"
    geo = folder + "geometry.in"
    num_points = 50
    geo = []
    with open(folder + "geometry.in", "r") as f:
        geo = f.readlines()
    for x in range(num_points):
        new_dir = folder + "relax_" + str(x)
        if x < 10:
            new_dir = folder + "relax_0" + str(x)
        angle = ((360 / num_points) * x) * (2 * math.pi) / 360
        x_coord = 4 * math.cos(angle)
        y_coord = 4 * math.sin(angle)
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        command = f'cp %scontrol.in %s' % (folder, new_dir)
        os.system(command)
        command = f'cp %ssubmit.sh %s' % (folder, new_dir)
        os.system(command)
        with open(new_dir + "/geometry.in", "w") as f:
            f.writelines(geo)
            f.write(f'atom\t\t%s\t\t%s\t\t0 I' % (x_coord, y_coord))


def move_FAI_geos():
    folder = "../../FHI-aims/French_NMR/FA/FAI/relax_many/"
    geo = folder + "geometries"
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dir.sort()
    os.mkdir(geo)
    i = 0
    for direct in all_dir:
        current = folder + direct
        if os.path.isdir(folder + direct):
            count = str(i)
            if i < 10:
                count = "0" + str(i)
            command = f'cp %s/geometry.in.next_step %s/geometry_%s.in' % (current, geo, count)
            os.system(command)
            i += 1


def constrain_base_relax(readfile, writefile):
    with open(readfile, "r") as f:
        lines = f.readlines()
    count = 0
    with open(writefile, "w") as f:
        for line in lines:
            if "lattice" in line:
                f.write(line)
            elif "Bi" in line:
                f.write(line)
                f.write("\t constrain_relaxation .true. \n")
            elif "Ag" in line:
                f.write(line)
                f.write("\t constrain_relaxation .true. \n")
            elif "I" in line:
                f.write(line)
                count += 1
                if count < 5:
                    f.write("\t constrain_relaxation x \n")
                    f.write("\t constrain_relaxation y \n")
                else:
                    f.write("\t constrain_relaxation .true. \n")
            elif "Pb" in line:
                f.write(line)
                f.write("\t constrain_relaxation .true. \n")
            elif "Tl" in line:
                f.write(line)
                f.write("\t constrain_relaxation .true. \n")
            elif "Br" in line:
                f.write(line)
                count += 1
                if count < 5:
                    f.write("\t constrain_relaxation x \n")
                    f.write("\t constrain_relaxation y \n")
                else:
                    f.write("\t constrain_relaxation .true. \n")
            elif "Cs" in line:
                f.write(line)
                f.write("\t constrain_relaxation x \n")
                f.write("\t constrain_relaxation y \n")
            else:
                print("Unrecognized species")


def make_Si_EV():
    folder = "../../FHI-aims/random/Si_EV/"
    base = folder + "base/"
    with open(base + "geometry.in", "r") as f:
        lines = f.readlines()
    for dist in ["5.25", "5.30", "5.35", "5.40", "5.45", "5.50", "5.55", "5.60", "5.65"]:
        temp = folder + dist + "/"
        if not os.path.exists(temp):
            os.mkdir(temp)
        command = "cp " + base + "control.in " + temp + "control.in"
        os.system(command)
        with open(temp + "geometry.in", "w") as f:
            f.write(f'lattice_vector %s 0 0\n' % dist)
            f.write(f'lattice_vector 0 %s 0\n' % dist)
            f.write(f'lattice_vector 0 0 %s\n' % dist)
            f.writelines(lines)


def exponential_func(x, a, b):
    return a * np.exp(b * x)


def ME431_lab1():
    pos = np.array([0, 60, 120, 180, 240, 300, 360])
    polished = np.array([75.53730117, 57.18318535, 42.52891284, 28.62567035, 20.26701245, 13.06427753, 11.18923853])
    painted = np.array([65.47410525, 47.08603857, 32.52502248, 21.47837276, 15.98674103, 10.52241235, 8.797268541])

    # Fit the exponential model
    params_pol, covariance_pol = curve_fit(exponential_func, pos, polished, p0=(1, -0.1))
    params_pain, covariance_pain = curve_fit(exponential_func, pos, painted, p0=(1, -0.1))

    # Extract parameters a and b
    a_pol, b_pol = params_pol
    a_pain, b_pain = params_pain

    # Standard deviations of the parameters (square root of the diagonal of covariance matrix)
    perr_pol = np.sqrt(np.diag(covariance_pol))
    a_pol_err, b_pol_err = perr_pol
    perr_pain = np.sqrt(np.diag(covariance_pain))
    a_pain_err, b_pain_err = perr_pain

    # Compute 95% confidence intervals (assuming a normal distribution, 1.96 * standard error gives approx. 95%)
    confidence_interval = 1.96

    # Predict y values for the fit
    y_fit_pol = exponential_func(pos, a_pol, b_pol)
    y_fit_pain = exponential_func(pos, a_pain, b_pain)

    # Upper and lower bounds of the 95% confidence interval for the y-values
    y_fit_upper_pol = exponential_func(pos, a_pol + confidence_interval * a_pol_err,
                                       b_pol + confidence_interval * b_pol_err)
    y_fit_lower_pol = exponential_func(pos, a_pol - confidence_interval * a_pol_err,
                                       b_pol - confidence_interval * b_pol_err)
    y_fit_upper_pain = exponential_func(pos, a_pain + confidence_interval * a_pain_err,
                                        b_pain + confidence_interval * b_pain_err)
    y_fit_lower_pain = exponential_func(pos, a_pain - confidence_interval * a_pain_err,
                                        b_pain - confidence_interval * b_pain_err)

    # Plot the original data, best fit, and the 95% confidence interval
    plt.scatter(pos, polished, label='Polished Rod')
    plt.plot(pos, y_fit_pol, label=f'Fit: y = {a_pol:.2f} * exp({b_pol:.2f} * x)', color='red')

    plt.scatter(pos, painted, label='Painted Rod')
    plt.plot(pos, y_fit_pain, label=f'Fit: y = {a_pain:.2f} * exp({b_pain:.2f} * x)', color='red')

    # Plot the confidence intervals
    plt.fill_between(pos, y_fit_lower_pol, y_fit_upper_pol, color='gray', alpha=0.3, label='_')
    plt.fill_between(pos, y_fit_lower_pain, y_fit_upper_pain, color='gray', alpha=0.3, label='_')

    # Add legend and labels
    plt.legend()
    plt.xlabel('Distance from heater (mm)')
    plt.ylabel('Temperature above ambient (C)')
    plt.title('Exponential Fit with 95% Confidence Interval')
    plt.show()

    print(f"a = {a_pol} ± {a_pol_err}")
    print(f"b = {b_pol} ± {b_pol_err}")

    print(f"a = {a_pain} ± {a_pain_err}")
    print(f"b = {b_pain} ± {b_pain_err}")