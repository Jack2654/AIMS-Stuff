import PlottingTools as pt
import Geometry
import KellerPBE as kp
import BasicGeo as bg
import BasicFunc as bf
import VectorToolkit as vt

debug = True

if not debug:
    file = "../../FHI-aims/KellerPBE/dissociation_curves/min_s_ts/"
    control = "../../FHI-aims/KellerPBE/control_files/control_pbe_ts.in"
    base_defaults = "../../FHI-aims/Repository/species_defaults/defaults_2020/tight/"
    min_defaults = "../../FHI-aims/Repository/species_defaults/min_s_defaults/"
    # kp.write_controls_to_dc(file, control, defaults)
if debug:
    file = "../../FHI-aims/KellerPBE/dissociation_curves/pbe_tight_ts/"
    kp.run_all_dc(file)