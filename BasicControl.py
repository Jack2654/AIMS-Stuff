# This file compiles basic functions for doing shit with control.in files
# In particular, the following functions exist here:
#
# -> write_control(geometry, base_control, defaults):
#                   given a geometry.in file, a template control.in, and a folder with species
#                   defaults, writes control.in file in same directory as geometry.in file
# -> write_all_controls(directory, base_control, defaults):
#                   using base_control as a template, writes a control file in each directory in the given directory
# -> species_default_path(element, defaults_folder):
#                   given an element and a path to the defaults, returns the path to the species default file
# -> list_of_defaults(filepath, defaul):
# ->                returns a list of paths to all required defaults files for the given folder
# ->

# imports:
import os
import BasicFunc as bf
import BasicGeo as bg
from PyAstronomy import pyasl


# functions:
def write_control(geometry, base_control, defaults):
    if not os.path.exists(geometry):
        print("bad path given to write_control " + geometry)
        return 0
    defaults = list_of_defaults(geometry, defaults)
    with open(base_control, "r") as f:
        cont = f.read()
    with open(geometry[:-11] + "control.in", "w") as f:
        f.writelines(cont)
        for default in defaults:
            with open(default, "r") as g:
                f.write(g.read())


def write_all_controls(directory, base_control, defaults):
    for direct in os.listdir(os.fsencode(directory)):
        dir_name = os.fsdecode(direct)
        geometry = directory + dir_name + "/geometry.in"
        write_control(geometry, base_control, defaults)


def species_default_path(element, defaults_folder):
    an = pyasl.AtomicNo()
    number = f"{an.getAtomicNo(element):02}"
    return defaults_folder + number + "_" + element + "_default"


def list_of_defaults(filepath, defaults_folder):
    defaults = []
    elements = bg.get_species_from_geo(filepath)
    for e in elements:
        defaults.append(species_default_path(e, defaults_folder))
    return defaults
