import os
import VectorToolkit as vt
import BasicGeo
import BasicFunc as bf
import BasicControl
from mins import Species
from molmass import Formula


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


def create_min_s_defaults(filepath):
    elem_dict = {'H': '1', 'He': '2', 'Li': '3', 'Be': '4', 'B': '5', 'C': '6', 'N': '7', 'O': '8',
            'F': '9', 'Ne': '10', 'Na': '11', 'Mg': '12', 'Al': '13', 'Si': '14', 'P': '15',
            'S': '16', 'Cl': '17', 'Ar': '18', 'K': '19', 'Ca': '20', 'Sc': '21', 'Ti': '22',
            'V': '23', 'Cr': '24', 'Mn': '25', 'Fe': '26', 'Co': '27', 'Ni': '28', 'Cu': '29',
            'Zn': '30', 'Ga': '31', 'Ge': '32', 'As': '33', 'Se': '34', 'Br': '35', 'Kr': '36',
            'Rb': '37', 'Sr': '38', 'Y': '39', 'Zr': '40', 'Nb': '41', 'Mo': '42', 'Tc': '43',
            'Ru': '44', 'Rh': '45', 'Pd': '46', 'Ag': '47', 'Cd': '48', 'In': '49', 'Sn': '50',
            'Sb': '51', 'Te': '52', 'I': '53', 'Xe': '54', 'Cs': '55', 'Ba': '56', 'La': '57',
            'Ce': '58', 'Pr': '59', 'Nd': '60', 'Pm': '61', 'Sm': '62', 'Eu': '63', 'Gd': '64',
            'Tb': '65', 'Dy': '66', 'Ho': '67', 'Er': '68', 'Tm': '69', 'Yb': '70', 'Lu': '71',
            'Hf': '72', 'Ta': '73', 'W': '74', 'Re': '75', 'Os': '76', 'Ir': '77', 'Pt': '78',
            'Au': '79', 'Hg': '80', 'Tl': '81', 'Pb': '82', 'Bi': '83', 'Po': '84', 'At': '85',
            'Rn': '86', 'Fr': '87', 'Ra': '88', 'Ac': '89', 'Th': '90', 'Pa': '91', 'U': '92',
            'Np': '93', 'Pu': '94', 'Am': '95', 'Cm': '96', 'Bk': '97', 'Cf': '98', 'Es': '99',
            'Fm': '100', 'Md': '101', 'No': '102'}
    for x in elem_dict.keys():
        species = Species(x)
        species.write_file(filepath)


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


def find_total_energy(aims):
    with open(aims, "r") as f:
        for line in f:
            if "| Total energy of" in line:
                return line.split()[11]
    return "ERROR: no total energy found in" + aims


def binding_energy(folder, output_file):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
    total_energies = {}
    binding_energies = {}
    monomers = 0
    for file in all_dir:
        if os.path.isdir(file):
            total_energies[file] = find_total_energy(folder + "/" + file + "/aims.out")
            if "monomer" in file:
                monomers += float(total_energies[file])
    for x in total_energies:
        if "monomer" not in x:
            binding_energies[x] = float(total_energies[x]) - monomers
    sorted_be = sorted(binding_energies)
    with open(output_file, "a") as f:
        f.write(folder[-2:] + "\n")
        for x in sorted_be:
            f.write(str(binding_energies[x]) + "\n")


def binding_energies(base, output):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        binding_energy(base + dir, output)


def generate_info(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    total_energies = {}
    for file in all_dir:
        if os.path.isdir(base + file):
            total_energies[file] = find_total_energy(base + "/" + file + "/aims.out")
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


def mass(element):
    return Formula(element).nominal_mass


def center_of_mass(base, target, even=False, a=False):
    tar = []

    if target:
        at = BasicGeo.atoms(base)
        for t in target:
            tar.append(at[t-1])
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
            cur_mass = mass(temp[4])
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


def write_controls_to_dc(base, control, defaults):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS" not in dir:
            BasicControl.write_all_controls(base + dir + "/", control, defaults)


def run_all_dc(base):
    all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(base))]
    all_dir.sort()
    for dir in all_dir:
        if ".DS" not in dir:
             bf.run_all(base + dir + "/")
