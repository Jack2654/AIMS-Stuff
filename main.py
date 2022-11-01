import math

import numpy as np
from math import pi
import Rotation
import Geometry
import m353

b = 170
a = 150

tilt = "../../FHI-aims/TlBi-Perovskites/out_of_plane/170/170_170/geometry.in"
output = "../../FHI-aims/TlBi-Perovskites/out_of_plane/170/170_" + str(a) + "/geometry.in"

temp = "./temp.in"

#Geometry.angleOOP(tilt, output, 1, 5, 8, b, a, True, True)
#Geometry.angleOOP(output, output, 1, 4, 17, b, a, True, False)
Geometry.angleOOP(output, temp, 17, 18, 13, b, a, False, True)
#Geometry.angleOOP(temp, output, 8, 9, 13, b, a, False, False)

theta = pi / 36
axis = np.array([0, 0, 1])
axis2 = np.array([0, 1, 0])
axis3 = np.array([1, 0, 0])

one = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/one.in"
two = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/two.in"
three = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/three.in"
four = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/four.in"
temp = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/temp.in"
temp2 = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/temp2.in"
built = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/built.in"
tempTl = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/tempTl.in"
tempBi = "../../FHI-aims/TlBi-Perovskites/setup_base/partials/tempBi.in"

# Rotation.recenter(one, one, 1, 0, 0, 0)

# Rotation.rotate(one, temp, theta, axis3)
# Ag
# Rotation.rotate(temp, temp2, theta, axis2)

#Rotation.rotate(three, temp, -theta, axis3)
# Bi
#Rotation.rotate(temp, temp2, -theta, axis2)

#Geometry.combine(tempTl, tempBi, built, 5, 2)
#Geometry.combine(built, tempTl, built, 10, 3)
#Geometry.combine(built, tempBi, built, 4, 3)
#Rotation.recenter(built, built, 1, 0, 0, 0)


fstart = "../../FHI-aims/TlBi-Perovskites/in_plane/"
fend = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/final_IP/180/180_190/band1001.out"

squirt = 1 / math.sqrt(2)
# Geometry.printBands(fstart, 1, .5)
# Geometry.printBands(fstart, 2, .5)
# Geometry.printBands(fstart, 3, squirt)
# Geometry.printBands(fstart, 4, squirt)


aang = 200
read = "../../FHI-aims/AgBi-Perovskites/ideal/tilting/final_IP/150/150_150/geometry.in"
write = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_vert/180/180_" + str(aang) + "/geometry.in"

# Geometry.addCs(read, read, 3)

# Geometry.delta_d_calc(read, 1, 2, 3, 4, 5, 6, 7, True)
