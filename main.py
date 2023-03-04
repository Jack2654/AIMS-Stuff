import math

from math import pi
import Rotation
import Geometry
import numpy as np
import matplotlib.pyplot as plt
import SpinTextures as st
import Methods2023 as m23
import egr244
import egr244 as dyna
import egr224 as mecha
import m566 as num
import PlottingTools as pt

b = 150
a = 170

# aang = 150
# newfile = "../../FHI-aims/TlBi-Perovskites/in_plane/"
# newfile += str(aang) + "/" + str(aang) + "_"
# for x in range(9):
#   fill = newfile + str(aang - 20 + x * 5) + "/geometry.in"
#   Geometry.centroid(fill, -10, 4, -10, 4, "Tl", False)

geo = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.next_step"
geo_two = "../../FHI-aims/ChiralSurfaceDefects/S_AEA/S_AEA_slab.in"
source = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceA2.in"
source2 = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceT2.in"
write = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.moved"
geo0 = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_OOP/180/180_200/band1001.out"
geo1 = "../../FHI-aims/AgBi-Perovskites/ideal/out/bro_idk/180/180_200/band1001.out"

# m23.manyLocAngles(source, write, 1, "s", True)
# print()
# m23.manyAbsAngles(source, write, 1, "s", True)
# m23.moveManyPoints(source, geo, 0, write, True)
# m23.moveManyPoints(source2, write, 0, write, False)
# m23.moveManyPoints(source2, write, 1, write, True)
# mecha.Delt_Y_Transform("D", 12, 6, 14)

geom_start_IP = "../../FHI-aims/AgBi-Perovskites/ideal/in_plane/final_IP/160/160_"
geom_start_OOP = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_OOP/160/160_"
geom_start_exp = "../../FHI-aims/AgBi-Perovskites/RFiles/bandsNEW/band100"
start = "../../FHI-aims/CationDoping/Final_Results/Doped_Relaxed/Zoom_new_control"
color_dict = {"Pb": "b", "I": "r", "C": "k", "H": "k", "N": "k", "Br": "k"}

#m23.remove_stars(start, start)
#m23.mullikenPlotter("../../FHI-aims/CationDoping/Final_Results/Doped_Experimental/Full_new_control/")

pt.mulliken_plot(start, start + "/out.png", "Pb", "I", 0.0, 2.1, 2.2, 0.0, color_dict, False)