import os
import math
import statistics
import VectorToolkit as vt
import BasicGeo
import BasicFunc as bf
import BasicControl
import BasicAimsOut as bao
import matplotlib.pyplot as plt
from scipy.optimize import Bounds
import time


def turn_all_to_geo(directory):
    files = []
    for file in os.listdir(os.fsencode(directory)):
        filename = os.fsdecode(file)
        files.append(filename)
    for f in files:
        read = directory + f
        write = "../../FHI-aims/KellerPBE/s66x10_in/" + f[:-4] + "/geometry.in"
        os.makedirs("../../FHI-aims/KellerPBE/s66x10_in/" + f[:-4])
        BasicGeo.xyz_to_geo(read, write)


def group_like_files(base):
    all_dir_b = os.listdir(os.fsencode(base))
    all_dir = []
    for x in all_dir_b:
        name = os.fsdecode(x)
        if "." not in name[0] and not len(name) == 10:
            all_dir.append(name)
    all_dir.sort()
    all_bases = set([x[:10] for x in all_dir])
    for b in all_bases:
        if not os.path.exists(base + b):
            os.makedirs(base + b)
    for dir in all_dir:
        os.rename(base + dir, base + dir[:10] + "/" + dir)


def binding_energy(folder, output_file):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    total_energies = {}
    binding_energies = {}
    monomers = 0
    for file in all_dir:
        if os.path.isdir(folder + "/" + file):
            total_energies[file] = bao.find_total_energy(folder + "/" + file + "/aims.out")
            if "A" in file or "B" in file or "monomer" in file:
                monomers += float(total_energies[file])
    for x in total_energies:
        if "monomer" not in x and "A" not in x and "B" not in x:
            if "ERR" in total_energies[x]:
                binding_energies[x] = -10000
            else:
                binding_energies[x] = float(total_energies[x]) - monomers
    sorted_be = sorted(binding_energies)
    with open(output_file, "a") as f:
        f.write(folder[-7:] + "\n")
        for x in sorted_be:
            f.write(str(binding_energies[x]) + "\n")


def binding_energies(base, output):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS_Store" not in dir and ".txt" not in dir:
            binding_energy(base + dir, output)


def generate_info(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    total_energies = {}
    for file in all_dir:
        if os.path.isdir(base + file):
            total_energies[file] = bao.find_total_energy(base + "/" + file + "/aims.out")
    with open(base + "info.txt", "w") as f:
        f.write("Total energies per file:\n")
        for x in total_energies.keys():
            f.write(x + ": " + total_energies[x] + "\n")
    print("Generated info file about total energies at " + base + "info.txt")


def formatS(name):
    temp = name[8:]
    if not temp[1].isdigit():
        temp = "0" + temp
    if "monomer" in temp:
        temp = temp[:3] + temp[10]
    else:
        if temp[6].isdigit():
            temp = temp[0:7]
        else:
            temp = temp[0:6] + "0"
    return temp


def better_geometries_folder(base, new):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    if not os.path.exists(new):
        os.makedirs(new)
    for direct in all_dir:
        command = "cp " + base + direct + "/geometry.in " + new + formatS(direct) + ".in"
        os.system(command)


def move_all_1_0_s(source, dest):
    all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(source))]
    all_geos.sort()
    if not os.path.exists(dest):
        os.makedirs(dest)
    for file in all_geos:
        if ("A" in file) or ("B" in file) or ("1.00" in file):
            command = "cp " + source + file + " " + dest + file
            os.system(command)


def organize_1_0_s(source):
    all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(source))]
    all_geos.sort()
    for file in all_geos:
        if not os.path.exists(source + file[:2]):
            os.makedirs(source + file[:2])
        command = "mv " + source + file + " " + source + file[:2] + "/" + file
        os.system(command)


def generate_structure(base, pointer):
    i = 0.8
    for x in range(21):
        if not f'{i:.2f}' == "1.00":
            factor = i - 1
            dist = [factor * y for y in pointer]
            BasicGeo.move(base + base[-3:-1] + "_B.in", base + "temp.in", dist)
            with open(base + base[-3:-1] + f'_{i:.2f}.in', "w") as f:
                f.writelines(BasicGeo.atoms(base + base[-3:-1] + "_A.in"))
                f.writelines(BasicGeo.atoms(base + "temp.in"))
        i += 0.02


def center_of_mass(base, target, even=False, a=False):
    tar = []

    if target:
        at = BasicGeo.atoms(base)
        for t in target:
            tar.append(at[t - 1])
    else:
        if a:
            at = BasicGeo.atoms(base[:-7] + "A.in")
            tar = at
        else:
            at = BasicGeo.atoms(base)
            temp = len(BasicGeo.atoms(base[:-7] + "A.in"))
            i = 0
            for a in at:
                if i >= temp:
                    tar.append(a)
                i += 1

    average = [0, 0, 0]
    total_mass = 0
    for atom in tar:
        temp = atom.split()
        if not even:
            cur_mass = bf.mass(temp[4])
        else:
            cur_mass = 1
        total_mass += cur_mass
        average[0] += cur_mass * float(temp[1])
        average[1] += cur_mass * float(temp[2])
        average[2] += cur_mass * float(temp[3])
    average = [x / total_mass for x in average]
    return average


def generate_many_structures():
    structures = []
    full = []
    base = "../../FHI-aims/KellerPBE/new_geos/"
    for line in read_info_about_s66():
        if line[3] == 'True':
            structures.append([line[0], line[1].split(","), line[2].split(","), True])
        else:
            structures.append([line[0], line[1].split(","), line[2].split(","), False])
    for s in structures:
        if s[1][0] == "a":
            full.append([f'{eval(s[0]):02d}', [], [], s[3]])
        else:
            full.append([f'{eval(s[0]):02d}', [eval(x) for x in s[1]], [eval(x) for x in s[2]], s[3]])

    for line in full:
        number = line[0]
        a = line[1]
        b = line[2]
        h_bond = line[3]
        x1 = center_of_mass(base + number + "/" + number + "_1.00.in", a, h_bond, True)
        x2 = center_of_mass(base + number + "/" + number + "_1.00.in", b, h_bond, False)
        pointer = vt.vecSub(x2, x1)
        generate_structure(base + number + "/", pointer)


def check_geometries():
    structures = []
    full = []
    for line in read_info_about_s66():
        if line[3] == 'True':
            structures.append([line[0], line[1].split(","), line[2].split(","), True])
        else:
            structures.append([line[0], line[1].split(","), line[2].split(","), False])
    for s in structures:
        if s[1][0] == "a":
            full.append([f'{eval(s[0]):02d}', [], [], s[3]])
        else:
            full.append([f'{eval(s[0]):02d}', [eval(x) for x in s[1]], [eval(x) for x in s[2]], s[3]])
    with open("../../FHI-aims/KellerPBE/info/real_IMD.txt", "w") as f:
        f.write("ideal:\n")
        f.write("2.00\n")
        f.write("1.50\n")
        f.write("1.25\n")
        f.write("1.10\n")
        f.write("1.05\n")
        f.write("0.95\n")
        f.write("0.90\n")
        f.write("0.80\n")
        f.write("0.70\n")
    for line in full:
        number = line[0]
        a = line[1]
        b = line[2]
        h_bond = line[3]

        file1 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_1.00.in"
        file2 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_2.00.in"
        file3 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_1.50.in"
        file4 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_1.25.in"
        file5 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_1.10.in"
        file6 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_1.05.in"
        file7 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_0.95.in"
        file8 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_0.90.in"
        file9 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_0.80.in"
        file10 = "../../FHI-aims/KellerPBE/s66x10_geos/xx_0.70.in"
        file1 = file1.replace("xx", number)
        file2 = file2.replace("xx", number)
        file3 = file3.replace("xx", number)
        file4 = file4.replace("xx", number)
        file5 = file5.replace("xx", number)
        file6 = file6.replace("xx", number)
        file7 = file7.replace("xx", number)
        file8 = file8.replace("xx", number)
        file9 = file9.replace("xx", number)
        file10 = file10.replace("xx", number)

        with open("../../FHI-aims/KellerPBE/info/real_IMD.txt", "a") as f:
            f.write(number + "\n")
            f.write(ratio(file1, a, file2, b, h_bond) + "\n")
            f.write(ratio(file1, a, file3, b, h_bond) + "\n")
            f.write(ratio(file1, a, file4, b, h_bond) + "\n")
            f.write(ratio(file1, a, file5, b, h_bond) + "\n")
            f.write(ratio(file1, a, file6, b, h_bond) + "\n")
            f.write(ratio(file1, a, file7, b, h_bond) + "\n")
            f.write(ratio(file1, a, file8, b, h_bond) + "\n")
            f.write(ratio(file1, a, file9, b, h_bond) + "\n")
            f.write(ratio(file1, a, file10, b, h_bond) + "\n")


def ratio(file1, a, file2, b, h_bond):
    one = vt.dist(center_of_mass(file1, a, h_bond, True), center_of_mass(file1, b, h_bond))
    two = vt.dist(center_of_mass(file2, a, h_bond, True), center_of_mass(file2, b, h_bond))
    return str(two / one)


def read_info_about_s66():
    with open("../../FHI-aims/KellerPBE/info/dissociation.dat", "r") as f:
        structures = []
        for line in f:
            structures.append(line.split())
    return structures


# took directory with all inputs written from generate_many_structures and set them up to be run
def setup_dissociation_curves(base):
    all_dirs = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dirs.sort()
    for dir in all_dirs:
        if dir == ".DS_Store":
            continue
        else:
            organize_dc(base + dir + "/")


def organize_dc(base):
    all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_geos.sort()
    for geo in all_geos:
        if ("temp" not in geo) and ("in" in geo):
            if not os.path.exists(base + geo[:-3]):
                os.makedirs(base + geo[:-3])
            command = "mv " + base + geo + " " + base + geo[:-3] + "/geometry.in"
            os.system(command)


# write control.in files to the directories containing 1320 models across 6 structures for the new,
# extended dissociation curves
def write_controls_to_dc(base, control, defaults, additional=""):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS" not in dir and "txt" not in dir:
            BasicControl.write_all_controls(base + dir + "/", control, defaults, additional)


# runs run_all in each directory found
def run_all_dc(base, ignore=False):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS" not in dir and ".txt" not in dir:
            bf.run_all(base + dir + "/", ignore)


def determine_all_minima_s66(filepath, option=1):
    with open(filepath, "r") as f:
        lines = []
        for ln in f:
            lines.append(ln)
    x = []
    for i in range(21):
        x.append(0.8 + i * 0.02)
    results = []

    for i in range(66):
        y = [float(x) for x in lines[1 + i * 22: (i + 1) * 22]]
        if option == 1:
            results.append(bf.minimize_func(bf.fit_poly(x, y, 4), bnds=Bounds(lb=0.8, ub=1.2)))
        if option == 2:
            results.append(min_lowest_quad(x, y))
    return results


def determine_all_minima_d442(filepath, option=1):
    with open(filepath, "r") as f:
        lines = []
        for ln in f:
            lines.append(ln)
    x = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.25, 1.50, 2.00]
    results = []
    for i in range(442):
        y = [float(x) for x in lines[1 + i * 11: (i + 1) * 11]]
        # , bnds=Bounds(lb=-0.5, ub=1.2)
        if option == 1:
            ret = bf.minimize_func(bf.fit_poly(x, y, 5))
            results.append(ret)
        elif option == 2:
            # print(i)
            ret = min_lowest_quad(x, y)
            results.append(ret)
        else:
            print("ERROR: bad option setting passed to determine_all_minima_d442")
    return results


def min_lowest_quad(x, y):
    res = bf.minimize_func(lowest_quad(x, y))
    return res


def lowest_quad(x, y):
    points = []
    i = 0
    for elem in y:
        points.append([x[i], elem])
        i += 1
    min = 10
    min_dex = -1
    for count, point in enumerate(points):
        if point[1] < min:
            min = point[1]
            min_dex = count
    if min_dex == 0:
        x_min = [points[0][0], points[1][0], points[2][0]]
        y_min = [points[0][1], points[1][1], points[2][1]]
    elif min_dex == len(points) - 1:
        x_min = [points[len(points) - 3][0], points[len(points) - 2][0], points[len(points) - 1][0]]
        y_min = [points[len(points) - 3][1], points[len(points) - 2][1], points[len(points) - 1][1]]
    else:
        x_min = [points[min_dex - 1][0], points[min_dex][0], points[min_dex + 1][0]]
        y_min = [points[min_dex - 1][1], points[min_dex][1], points[min_dex + 1][1]]

    res = bf.fit_poly(x_min, y_min, 2)
    return res


def mean_absolute_error(file1, file2):
    r1 = determine_all_minima_s66(file1)
    r2 = determine_all_minima_s66(file2)
    error = []
    for i in range(len(r1)):
        error.append(abs(r1[i] - r2[i]))
    return sum(error) / len(error)


def plot_dissociation_curve_s66(files, structure):
    series = [f + "binding_energies.txt" for f in files]
    x = []
    for i in range(21):
        x.append(0.8 + i * 0.02)
    colors = ["red", "orange", "blue", "black", "purple"]
    i = 0
    for s in series:
        with open(s, "r") as f:
            temp = f.readlines()
        vals = [float(x) for x in temp[structure * 22 - 21: structure * 22]]
        plt.plot(x, vals, color=colors[i])
        plt.axvline(x=bf.minimize_func(bf.fit_poly(x, vals, 4), bnds=Bounds(lb=0.8, ub=1.2)), color=colors[i],
                    linestyle="--", label="_nolegend_")
        if False:
            print(bf.minimize_func(bf.fit_poly(x, vals, 4), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 5), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 6), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 7), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 8), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 9), bnds=Bounds(lb=0.8, ub=1.2)))
            print(bf.minimize_func(bf.fit_poly(x, vals, 10), bnds=Bounds(lb=0.8, ub=1.2)))
        i += 1
    plt.legend(["PBE tight", "PBE tight+ts", "min+s", "blyp_min"], loc="upper right")
    plt.title("Structure " + str(structure) + " Dissociation Curves")
    plt.show()


def plot_dissociation_curve_d442(files, structure):
    series = [f + "binding_energies.txt" for f in files]
    x = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.25, 1.50, 2.00]
    # colors = ["red", "orange", "blue", "black", "purple"]
    i = 0
    for s in series:
        with open(s, "r") as f:
            temp = f.readlines()
        vals = [float(x) for x in temp[structure * 11 - 10: structure * 11]]
        plt.plot(x, vals)
        # plt.axvline(x=min_lowest_quad(x, vals), linestyle="--", label="_nolegend_")
        i += 1
    legends = [f.split("/")[6] for f in files]
    plt.legend(legends)
    plt.title("D442 Structure " + str(structure) + " Dissociation Curves")
    plt.show()


def eq_form(coeff):
    res = ""
    i = 0
    for c in coeff:
        res += str(c) + "x^" + str(i) + "+"
        i += 1
    return res[:-1]


def s66x21_run(poss, base, control_file, defaults_folder, ignore=False, write_control=True):
    times = []
    for file in poss:
        print("Beginning calculations for settings:" + str(file))
        temp_time = time.time()
        temp_sr = str(file[0])
        temp_d = str(file[1])
        cur_folder = base + temp_sr + "_" + temp_d + "/"
        commands = "vdw_damping_sr " + temp_sr + "\n vdw_damping_d " + temp_d + "\n"
        be_file = cur_folder + "binding_energies.txt"

        if write_control:
            write_controls_to_dc(cur_folder, control_file, defaults_folder, commands)
        run_all_dc(cur_folder, ignore=ignore)
        binding_energies(cur_folder, be_file)

        re_time = time.time() - temp_time
        times.append(re_time)
        print(str(file) + " took " + str(re_time) + " seconds")

    for t in times:
        print("Operation took " + str(t) + " seconds")

    compare_to = "../../FHI-aims/KellerPBE/S66/dissociation_curves/pbe_tight_ts/binding_energies.txt"
    for file in poss:
        temp_sr = str(file[0])
        temp_d = str(file[1])
        cur_folder = base + temp_sr + "_" + temp_d + "/"
        be_file = cur_folder + "binding_energies.txt"
        print(str(file) + ": " + str(mean_absolute_error(be_file, compare_to)))

    for file in poss:
        temp_sr = str(file[0])
        temp_d = str(file[1])
        cur_folder = base + temp_sr + "_" + temp_d + "/"
        be_file = cur_folder + "binding_energies.txt"
        res = determine_all_minima_s66(be_file, 2)
        print()
        for r in res:
            print(r)


def organize_d442_1(folder):
    all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    for geo in all_geos:
        if not os.path.exists(folder + geo[:7]):
            os.makedirs(folder + geo[:7])
        command = "mv " + folder + geo + " " + folder + geo[:7] + "/."
        os.system(command)


def organize_d442_2(folder):
    all_dirs = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dirs.sort()
    for dir in all_dirs:
        all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder + dir + "/"))]
        all_geos.sort()
        for geo in all_geos:
            input = folder + dir + "/" + geo
            output = folder + dir + "/" + geo[:-3] + "in"
            BasicGeo.xyz_to_geo(input, output)
            command = "rm " + input
            os.system(command)


def organize_d442_3(folder):
    all_dirs = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    all_dirs.sort()
    for dir in all_dirs:
        all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder + dir + "/"))]
        all_geos.sort()
        for geo in all_geos:
            command = "mv " + folder + dir + "/" + geo + " " + folder + dir + "/" + geo[-6:-3] + "/geometry.in"
            os.makedirs(folder + dir + "/" + geo[-6:-3])
            os.system(command)


def d442x10_run(poss, base, commands, control_file, default_files, ignore=False, write_control=True):
    times = []

    i = -1
    for file in poss:
        i += 1
        temp_time = time.time()
        cur_folder = base + file
        be_file = cur_folder + "binding_energies.txt"
        if os.path.exists(be_file) and ignore is False and False:
            re_time = time.time() - temp_time
            times.append(re_time)
            print(str(file) + "\twas skipped and took " + str(re_time) + " seconds")
            continue

        print("Beginning calculations for: " + str(file))
        if write_control:
            write_controls_to_dc(cur_folder, control_file, default_files[i], commands[i])
        run_all_dc(cur_folder, ignore=ignore)
        binding_energies(cur_folder, be_file)

        re_time = time.time() - temp_time
        times.append(re_time)
        print(str(file) + " took " + str(re_time) + " seconds")

    for t in times:
        print("Operation took " + str(t) + " seconds")

    for file in poss:
        be_file = base + file + "binding_energies.txt"
        res = determine_all_minima_d442(be_file, option=2)
        print()
        print(file)
        for r in res:
            print(r)


def total_energy(folder, output_file):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    total_energies = {}
    for file in all_dir:
        if os.path.isdir(folder + "/" + file):
            total_energies[file] = bao.find_total_energy(folder + "/" + file + "/aims.out")
    sorted_te = sorted(total_energies)
    with open(output_file, "a") as f:
        f.write(folder[-7:] + "\n")
        for x in sorted_te:
            f.write(str(total_energies[x]) + "\n")


def all_energies(base, output):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS_Store" not in dir and ".txt" not in dir:
            total_energy(base + dir, output)


def split_monomers(input, output_base):
    with open(input, "r") as f:
        lines = f.readlines()
    monomers = lines[1].split()[3:5]
    monomer_a = [int(x) for x in monomers[0].split("=")[1].split("-")]
    monomer_b = [int(x) for x in monomers[1].split("=")[1].split("-")]

    if len(monomer_a) == 1:
        monomer_a = [monomer_a[0], monomer_a[0]]
    if len(monomer_b) == 1:
        monomer_b = [monomer_b[0], monomer_b[0]]

    mon_a = []
    for count, line in enumerate(lines):
        if monomer_a[0] < count < monomer_a[1] + 2:
            temp = line.strip().split()
            temp.append(temp[0])
            temp[0] = "atom"
            mon_a.append("   ".join(temp) + "\n")

    mon_b = []
    for count, line in enumerate(lines):
        if monomer_b[0] < count < monomer_b[1] + 2:
            temp = line.strip().split()
            temp.append(temp[0])
            temp[0] = "atom"
            mon_b.append("   ".join(temp) + "\n")

    if len(mon_a) + len(mon_b) + 2 != len(lines):
        print("Monomer Split Error: " + input)
        return 0

    with open(output_base + "A/geometry.in", "w") as f:
        f.writelines(mon_a)
    with open(output_base + "B/geometry.in", "w") as f:
        f.writelines(mon_b)


def all_d442_monomers(base, output):
    all_geos = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_geos.sort()
    geos = []
    for geo in all_geos:
        geos.append(geo[:-8])
    geos = set(geos)
    for geo in geos:
        input = base + geo + "_100.xyz"
        output_base = output + geo + "/00"
        os.mkdir(output_base + "A")
        os.mkdir(output_base + "B")
        split_monomers(input, output_base)


def plot_s66_histogram(series, comparison):
    colors = ["green", "orange", "blue"]
    base_be = comparison + "binding_energies.txt"
    base_vals = determine_all_minima_s66(base_be)
    i = 0
    for file in series:
        comp_be = file + "binding_energies.txt"
        comp_vals = determine_all_minima_s66(comp_be)
        all_errors = []
        for count, val in enumerate(comp_vals):
            error = base_vals[count] - val
            all_errors.append(error)
        plt.hist(all_errors, bins=30, alpha=0.4, color=colors[i])
        plt.axvline(x=[statistics.median(all_errors)], label="_nolegend_", color=colors[i])
        # plt.axvline(x=[statistics.mean(all_errors)], label="_nolegend_", color=colors[i], linestyle="--")
        i += 1
    plt.axvline(x=[0.0], color='k', label="_nolegend_")
    plt.xlabel("Error (% of IMD)")
    plt.ylabel("Count")
    plt.title("S66x21 Structures")
    plt.legend(["min+s", "min+s + ts", "min+s sr=1.00, d=70"])
    plt.show()


def plot_d442_histogram(series, comparison):
    colors = ["green", "orange", "blue", "red"]
    base_be = comparison + "binding_energies.txt"
    base_vals = determine_all_minima_d442(base_be, option=2)
    i = 0
    for file in series:
        comp_be = file + "binding_energies.txt"
        comp_vals = determine_all_minima_d442(comp_be, option=2)
        all_errors = []
        for count, val in enumerate(comp_vals):
            error = base_vals[count] - val
            all_errors.append(error)
        errors = []
        for err in all_errors:
            if -0.5 < err < 0.5:
                errors.append(err)
        plt.hist(errors, bins=30, alpha=0.5, color=colors[i+2])
        plt.axvline(x=[statistics.median(all_errors)], label="_nolegend_", color=colors[i+2])
        # plt.axvline(x=[statistics.mean(all_errors)], label="_nolegend_", color=colors[i], linestyle="--")
        i += 1
    plt.axvline(x=[0.0], color='k', label="_nolegend_")
    plt.xlabel("Error (% of IMD)")
    plt.ylabel("Count")
    plt.title("D442x10 Structures")
    plt.legend([x.split("/")[6] for x in series]) # was changed by accident, check if still works
    plt.show()
