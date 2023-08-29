from ase.calculators.aims import Aims
import os
from ase.calculators.mixing import SumCalculator
from ase.io import read
from mins import SmallBasisCorrection

mpirun = "/opt/homebrew/bin/orterun -N 8 "
aims_exec = "/Users/jackmorgenstein/FHI-aims/Repository/builds/build_06_09_23/aims.230612.mpi.x"
path = " > aims.out"
command = mpirun + aims_exec + path
species_tight = "/Users/jackmorgenstein/FHI-aims/Repository/species_defaults/defaults_2020/tight/"
species_min = "/Users/jackmorgenstein/FHI-aims/Repository/species_defaults/min_s_defaults/"

os.system("ulimit -s hard")
os.system("export TMPDIR=/tmp")
os.system("export OMP_NUM_THREADS=1")
os.environ['ASE_AIMS_COMMAND'] = command

calc_fhiaims_mins = Aims(xc="pbe", vdw_correction_hirshfeld=".true.",
                         relativistic="atomic_zora scalar", spin="none", check_stacksize=".false.")
calc_mins_correction = SmallBasisCorrection()
structure = read("../../FHI-aims/KellerPBE/ASE_test/test_1/geometry.in")
structure.calc = calc_fhiaims_mins


print("tight species:")
calc_fhiaims_mins.set(species_dir=species_tight)
print(structure.get_potential_energy())

print("min species:")
calc_fhiaims_mins.set(species_dir=species_min)
print(structure.get_potential_energy())

# calc = SumCalculator([calc_fhiaims_mins, calc_mins_correction])
structure = read("../../FHI-aims/KellerPBE/ASE_test/test_1/geometry.in")
# structure.calc = calc
# structure.get_potential_energy()

calc_mins_correction.calculate(atoms=structure)
print(calc_mins_correction.get_potential_energy())


# calc_mins_correction.calculate(atoms=structure)
# structure.calc = calc_fhiaims_mins
# print(structure.get_potential_energy())

