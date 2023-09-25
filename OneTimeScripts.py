# file for 1 time use scripts

# imports
import BasicGeo as bg
import os

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
                if "  species  " in line:
                    f.write(line)
                    f.write("  species_default_type min+s\n")
                else:
                    f.write(line)
