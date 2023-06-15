import PlottingTools as pt
import Geometry
import ASE_tools as ASE_t
import Basic
import VectorToolkit as vt

debug = False

if not debug:
    ASE_t.generate_many_structures()

if debug:
    ASE_t.check_geometries()