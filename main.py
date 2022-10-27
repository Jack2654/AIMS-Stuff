import math

import numpy as np
from math import pi
import Rotation
import Geometry
import m353


theta = pi / 12
axis = np.array([0, 0, 1])
axis2 = np.array([0, 1, 0])

one = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/one.in"
two = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/two.in"
three = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/three.in"
four = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/four.in"
temp = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/temp.in"
temp2 = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/temp2.in"
built = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/built.in"

#Rotation.rotate(one, temp, theta, axis)
#Rotation.rotate(two, temp2, -theta, axis)
#Geometry.combine(temp2, temp, built, 3, 4)
#Rotation.rotate(three, temp, -theta, axis)
#Geometry.combine(built, temp, built, 10, 2)
#Rotation.rotate(four, temp, theta, axis)
#Geometry.combine(built, temp, built, 15, 2)
#Rotation.recenter(built, built, 8, 0, 0, 0)


fileR = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_150/geometry.in"
fileW = "../../FHI-aims/AgBi-Perovskites/ideal/out/150c/150_160/geometry.in"
b = 150
a = 160

#Geometry.angleOOP(fileR, fileW, 1, 2, 12, b, a, True, True)
#Geometry.angleOOP(fileW, fileW, 1, 3, 7, b, a, True, False)
#Geometry.angleOOP(fileW, fileW, 6, 7, 15, b, a, False, True)
#Geometry.angleOOP(fileW, fileW, 9, 10, 14, b, a, False, False)

fstart = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_vert/"

squirt = 1/math.sqrt(2)
#Geometry.printBands(fstart, 3, squirt)
#Geometry.printBands(fstart, 4, squirt)

#Geometry.angleInfo("../../FHI-aims/AgBi-Perovskites/ideal/tilting/180/bands/180_180/geometry.in", 3, 25, 1, True)

f1 = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/final_IP/180/180_180/geometry.in"
f2 = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/final_IP/180/180_185/geometry.in"
