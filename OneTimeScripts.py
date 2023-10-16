# file for 1 time use scripts

# imports
import BasicGeo as bg
import os
import matplotlib.pyplot as plt

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