import math

import numpy as np
from math import pi
import Rotation
import Geometry
import m353

#m353.manyEuler()

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

filea = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_130/band1001.out"
fileb = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_135/band1001.out"
filec = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_140/band1001.out"
filed = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_145/band1001.out"
filee = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_150/band1001.out"
filef = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_155/band1001.out"
fileg = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_160/band1001.out"
fileh = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_165/band1001.out"
filei = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/150/bands2/150_170/band1001.out"

a = 160
fileI = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/160/bands2/160_"
fileO = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_"
for i in range(9):
    fileplus = fileI + str((a - 20)+(5*i))
    fileplus += "/band1001.out"
    Geometry.bandProcess(fileplus, True)
for i in range(9):
    fileplus = fileI + str((a - 20)+(5*i))
    fileplus += "/band1001.out"
    Geometry.bandProcess(fileplus, False)
for i in range(9):
    fileplus = fileI + str((a - 20)+(5*i))
    fileplus += "/band1001.out"
    Geometry.bandGap(fileplus)

# random comments
#Geometry.bandProcess(filea)
#Geometry.bandProcess(fileb)
#Geometry.bandProcess(filec)
#Geometry.bandProcess(filed)
#Geometry.bandProcess(filee)
#Geometry.bandProcess(filef)
#Geometry.bandProcess(fileg)
#Geometry.bandProcess(fileh)
#Geometry.bandProcess(filei)

#Geometry.bandGap(filea)
#Geometry.bandGap(fileb)
#Geometry.bandGap(filec)
#Geometry.bandGap(filed)
#Geometry.bandGap(filee)
#Geometry.bandGap(filef)
#Geometry.bandGap(fileg)
#Geometry.bandGap(fileh)
#Geometry.bandGap(filei)