def set_header(basis, elem):
    not_benchmarked = ["He",
                       "Li", "Be", "Ne",
                       "Na", "Mg", "Ar"]
    header = "################################################################################\n"
    header += "#\n"
    header += "#  FHI-aims code project\n"
    header += "#\n"
    header += f'#  Suggested \"{basis}\" defaults for {elem}\n'
    header += "#\n"
    header += "#  The NAO-J basis sets are intended for very tightly converged calculations,\n"
    header += "#  especially around the nucleus, to compute numerically precise spin-spin coupling\n"
    header += "#  constants (\"J-couplings\") for NMR.\n"
    header += "#\n"
    header += "#  Definitions and benchmarks are provided in\n"
    header += "#\n"
    header += "#  Laasner et al., Electronic Structure 6, 027002 (2024). Please see Figure 8 and\n"
    header += "#  associated discussion in that paper for details.\n"
    header += "#\n"
    if elem in not_benchmarked:
        header += f'#  The basis set for {elem} was NOT benchmarked in that paper and was constructed later\n'
        header += "#  by Jack Morgenstein and Volker Blum, applying the same construction principle.\n"
        header += "#\n"
    header += "################################################################################\n"
    header += f'    species        {elem}\n'
    header += "#     set basis set reference for FHI-aims output\n"
    header += "    cite_reference NAO-J-n-2024\n"
    header += "#     global species definitions\n"
    return header


def write_NAO_J():
    vcc_sets = ["NAO-VCC-2Z/", "NAO-VCC-3Z/", "NAO-VCC-4Z/", "NAO-VCC-5Z/"]
    j_sets = ["NAO-J-2", "NAO-J-3", "NAO-J-4", "NAO-J-5"]
    j_sets = ["NAO-J-p2", "NAO-J-p3", "NAO-J-p4", "NAO-J-p5"]
    species_arr = ['H', 'He',
                   'Li', 'Be', 'B',  'C',  'N', 'O', 'F',  'Ne',
                   'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

    directory = "/Users/jackmorgenstein/FHI-aims/Repositories/FHIaims_NAO_J/species_defaults/defaults_next/nmr/NAO-J-n/"
    # directory = "/Users/jackmorgenstein/FHI-aims/Repositories/FHIaims_NAO_J/species_defaults/defaults_next/nmr/NAO-J-pn/"
    for count_basis, basis_type in enumerate(vcc_sets):
        direct = f'/Users/jackmorgenstein/FHI-aims/Repositories/FHIaims_NAO_J/species_defaults/NAO-VCC-nZ/{basis_type}'
        for count_elem, elem in enumerate(species_arr):
            cutoff = 1000
            if elem == "H":
                cutoff = 100

            with open(direct + f'{str(count_elem+1).rjust(2,"0")}_{elem}_default', "r") as f:
                vcc_lines = f.readlines()

            with open(directory + f'../pcJ-4_temp/{str(count_elem+1).rjust(2,"0")}_{elem}', "r") as f:
                pc_lines = f.readlines()

            header = set_header(j_sets[count_basis], elem)
            with open(directory + j_sets[count_basis] + f'/{str(count_elem+1).rjust(2,"0")}_{elem}_default', "w") as f:
                f.write(header)
                started = False
                for line in vcc_lines:
                    if "global species" in line:
                        started = True
                    elif started:
                        if "l_hartree" in line:
                            f.write("    l_hartree           8\n")
                        elif "radial_mult" in line:
                            f.write("    radial_multiplier   8\n")
                        elif "basis_dep_cutoff" in line:
                            f.write(line)
                            f.write("    logarithmic         0.00001 100 1.0123\n")
                        else:
                            f.write(line)
                searching = False
                for line in pc_lines:
                    if "gaussian" in line:
                        searching = False
                        temp = line.split()
                        if temp[1] == "0" or temp[1] == "1":  # !!!!!!! change this line to ensure only s-type orbitals
                            if temp[2] != "1":
                                searching = True
                            else:
                                if temp[1] == "0" and float(temp[3]) > cutoff:
                                    f.write(f' gaussian 0 1  {temp[3]}\n')
                                elif temp[1] == "1" and float(temp[3]) > 10:
                                    f.write(f' gaussian 1 1  {temp[3]}\n')
                    elif searching:
                        for x in line.split():
                            if float(x) > cutoff:
                                f.write(f' gaussian 0 1  {x}\n')


def write_pc_J(num):
    species_arr = ['H', 'He',
                   'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                   'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

    pcJ_direct = f'/Users/jackmorgenstein/FHI-aims/Repositories/FHIaims_pcJ/species_defaults/defaults_next/nmr/tmp/new_{num}/'
    vcc_direct = "/Users/jackmorgenstein/FHI-aims/Repositories/FHIaims_pcJ/species_defaults/NAO-VCC-nZ/NAO-VCC-5Z/"
    for count_elem, elem in enumerate(species_arr):
        with open(vcc_direct + f'{str(count_elem + 1).rjust(2, "0")}_{elem}_default', "r") as f:
            vcc_lines = f.readlines()

        with open(pcJ_direct + f'../old_{num}/{str(count_elem + 1).rjust(2, "0")}_{elem}', "r") as f:
            pc_lines = f.readlines()

        BSE_header = "################################################################################\n" + \
                     "# Basis Set Exchange\n" + \
                     "# Version 0.10\n" + \
                     "# https://www.basissetexchange.org\n" + \
                     "################################################################################\n" + \
                     f'# Basis set: pcJ-{num} for {elem}\n' + \
                     f'# Description: Contracted version of the pcJ-{num} basis\n' + \
                     "# Role: orbital\n" + \
                     "# Version: 1  (Data from Frank Jensen)\n" + \
                     "################################################################################\n" + \
                     "# Basis set completed in FHI-aims format by Jack Morgenstein and Volker Blum\n" + \
                     "################################################################################\n" + \
                     f'    species        {elem}\n'

        with open(pcJ_direct + f'/{str(count_elem + 1).rjust(2, "0")}_{elem}_default', "w") as f:
            f.write(BSE_header)
            started = False
            for line in vcc_lines:
                if "global species" in line:
                    started = True
                elif started:
                    if "l_hartree" in line:
                        f.write("    l_hartree           8\n")
                    elif "radial_mult" in line:
                        f.write("    radial_multiplier   8\n")
                    elif "basis_dep_cutoff" in line:
                        f.write(line)
                        f.write("    logarithmic         0.00001 100 1.0123\n")
                    elif "basis_acc" in line:
                        started = False
                    else:
                        f.write(line)

            searching = False
            for line in pc_lines:
                if "pcJ" in line:
                    continue
                elif "gaussian" in line:
                    searching = False
                    temp = line.split()
                    if temp[2] != "1":
                        searching = True
                    f.write(line)
                elif searching:
                    f.write(f'      {line}')
                elif "pure_gauss" in line:
                    f.write(line + "\n")
                else:
                    f.write(line)


def split_elem(input_file, output_dir):
    with open(input_file, "r") as f:
        lines = f.readlines()
    elem = [[] for x in range(19)]
    index = 0
    for line in lines:
        if "default minimal basis" in line:
            index += 1
        elem[index].append(line)
    species_arr = ['H', 'He',
                   'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                   'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
    for count, el in enumerate(elem[1:]):
        with open(f'{output_dir}{str(count+1).rjust(2,"0")}_{species_arr[count]}', "w") as f:
            f.writelines([x.strip() + "\n" for x in el])


for i in range(5):
    number = str(i)
    input_file = f'../../FHI-aims/Repositories/FHIaims_pcJ/species_defaults/defaults_next/nmr/tmp/input_{number}'
    output_dir = f'../../FHI-aims/Repositories/FHIaims_pcJ/species_defaults/defaults_next/nmr/tmp/old_{number}/'
    split_elem(input_file, output_dir)
    write_pc_J(number)
