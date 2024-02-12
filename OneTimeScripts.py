# file for 1 time use scripts

# imports
import BasicGeo as bg
import os
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
import math


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


def visualize_spin_texture_directions():
    file = "../../FHI-aims/Double_Perovskites/sample_spin_text.in"
    points = []
    with open(file, "r") as f:
        for line in f:
            temp = line.split()
            points.append([float(x) for x in temp[2:8]])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for point in points:
        ax.plot([point[0], point[3]],
                [point[1], point[4]],
                [point[2], point[5]], 'k')
    plt.show()


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

