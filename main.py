import math

from math import pi
import Rotation
import Geometry
import numpy as np
import matplotlib.pyplot as plt
import SpinTextures as st
import Methods2023 as m23
import egr244

b = 150
a = 170

# aang = 150
# newfile = "../../FHI-aims/TlBi-Perovskites/in_plane/"
# newfile += str(aang) + "/" + str(aang) + "_"
# for x in range(9):
#   fill = newfile + str(aang - 20 + x * 5) + "/geometry.in"
#   Geometry.centroid(fill, -10, 4, -10, 4, "Tl", False)

geo = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.next_step"
geo2 = "../../FHI-aims/ChiralSurfaceDefects/S_AEA/S_AEA_slab.in"
source = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceA2.in"
source2 = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceT2.in"
write = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.moved"
m23.manyLocAngles(source, write, 1, "s", True)
print()
m23.manyAbsAngles(source, write, 1, "s", True)
#m23.moveManyPoints(source, geo, 0, write, True)
#m23.moveManyPoints(source2, write, 0, write, False)
#m23.moveManyPoints(source2, write, 1, write, True)

# data = egr244.hw2(0.1)
# print(data)
