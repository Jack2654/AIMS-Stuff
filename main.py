import math

import numpy as np
from math import pi
import Rotation
import Geometry

theta = pi / 12
axis = np.array([1, 0, 0])
axis2 = np.array([0, 1, 0])

fileR = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_150/geometry.in"
fileW = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_160/geometry.in"

b = 150
a = 160

#Geometry.angleOOP(fileR, fileW, 1, 2, 12, b, a, True, True)
#Geometry.angleOOP(fileW, fileW, 1, 3, 7, b, a, True, False)
#Geometry.angleOOP(fileW, fileW, 6, 7, 15, b, a, False, True)
#Geometry.angleOOP(fileW, fileW, 9, 10, 14, b, a, False, False)

filea = "../../FHI-aims/AgBi-Perovskites/ideal/out/160/160_180/band1003.out"
fileb = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_165/band1003.out"
filec = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_160/band1003.out"
filed = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_155/band1003.out"
filee = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_150/band1003.out"
filef = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_145/band1003.out"
fileg = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_140/band1003.out"
fileh = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_135/band1003.out"
filei = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_130/band1003.out"

# Geometry.bandProcess(filea)
# Geometry.bandProcess(fileb)
Geometry.bandProcess(filec)
Geometry.bandProcess(filed)
# Geometry.bandProcess(filee)
# Geometry.bandProcess(filef)
# Geometry.bandProcess(fileg)
# Geometry.bandProcess(fileh)
#Geometry.bandProcess(filei)