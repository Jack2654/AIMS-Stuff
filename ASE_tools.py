import os
from PyAstronomy import pyasl
import time
from mins import Species


def xyz_to_geo_s66x10(filepath, writepath):
    with open(filepath, "r") as f:
        at = []
        i = 0
        for ln in f:
            if i > 1:
                at.append(ln)
            i += 1
    with open(writepath, "w+") as f:
        for atom in at:
            sp = atom.split()
            temp = "atom " + sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[0] + "\n"
            f.write(temp)
    if not len(open(filepath, 'r').readlines()) == (len(open(writepath, 'r').readlines()) + 2):
        print("ERROR: error in translating .xyz to .in for " + filepath)


def turn_all_to_geo(directory):
    files = []
    for file in os.listdir(os.fsencode(directory)):
        filename = os.fsdecode(file)
        files.append(filename)
    for f in files:
        read = directory + f
        write = "../../FHI-aims/KellerPBE/s66x10_in/" + f[:-4] + "/geometry.in"
        os.makedirs("../../FHI-aims/KellerPBE/s66x10_in/" + f[:-4])
        xyz_to_geo_s66x10(read, write)


def get_species_from_geo(filepath):
    with open(filepath, "r") as f:
        species = []
        for ln in f:
            species.append(ln.split()[4])
    return set(species)


def species_default_path(element, defaults_folder):
    an = pyasl.AtomicNo()
    number = f"{an.getAtomicNo(element):02}"
    return defaults_folder + number + "_" + element + "_default"


def list_of_defaults(filepath, defaul):
    defaults = []
    elements = get_species_from_geo(filepath)
    for e in elements:
        defaults.append(species_default_path(e, defaul))
    return defaults


def write_control(geometry, write, defaults):
    defaults = list_of_defaults(geometry, defaults)
    with open(write, "w") as f:
        f.write("xc \t pbe \n relativistic \t atomic_zora scalar \n spin \t none \n check_stacksize \t .false. \n")
        for default in defaults:
            with open(default, "r") as g:
                f.write(g.read())


def write_all_controls(directory, defaults):
    for direct in os.listdir(os.fsencode(directory)):
        dir_name = os.fsdecode(direct)
        geometry = directory + dir_name + "/geometry.in"
        control = directory + dir_name + "/control.in"
        write_control(geometry, control, defaults)


def run_AIMS(directory):
    temp = os.getcwd()
    os.chdir(directory)
    os.system("ulimit -s hard")
    os.system("export TMPDIR=/tmp")
    os.system("export OMP_NUM_THREADS=1")
    os.system("/opt/homebrew/bin/orterun -N 8 /Users/jackmorgenstein/FHI-aims/Repository/builds/build_06_09_23/aims.230612.mpi.x > aims.out 2>&1")
    os.chdir(temp)


def run_all(base, ignore=False):
    all_dir_b = os.listdir(os.fsencode(base))
    all_dir = []
    for x in all_dir_b:
        all_dir.append(os.fsdecode(x))
    all_dir.sort()
    for direct in all_dir:
        if os.path.isdir(base + direct):
            dir_name = direct
            directory = base + dir_name
            if not os.path.exists(directory + "/aims.out") or ignore:
                current = time.time()
                run_AIMS(directory)
                print("Completed calculation for " + dir_name + " in " + str(time.time() - current) + " seconds")


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
