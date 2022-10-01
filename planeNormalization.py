import math
from sympy import *
import os.path


# Basic cross product formula, a and b should be three element tuples
def crossProduct(a, b):
    return a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]


# Basic 3 element dot product. a and b should be three element tuples
def dotProduct(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# This function scales the output values f1, f2, f3 (contained in tuple r) such that they output a vector within the
# Brillouin zone. NOTE: this does not actually compute the brillouin zone, it merely uses an approximation by
# checking in 26 distinct directions
def scale(t, r, b1, b2, b3):
    print("b1: ", b1, ", b2: ", b2, ", b3: ", b3)
    s = 1
    n = (-1, 0, 1)
    for x in range(3):
        for y in range(3):
            for z in range(3):

                # This if check is necessary so the component in the direction (0,0,0) is not calculated
                if x == 1 and y == 1 and z == 1:
                    continue
                dir = (n[x] * b1[0] + n[y] * b2[0] + n[z] * b3[0],
                       n[x] * b1[1] + n[y] * b2[1] + n[z] * b3[1],
                       n[x] * b1[2] + n[y] * b2[2] + n[z] * b3[2])
                comp = component(dir, t)
                if comp > 0.5:
                    if s > (1 / (2 * comp)):
                        # print("max dir: ", dir)
                        s = 1 / (2 * comp)
    # print("s: ", s)
    return tuple(s * x for x in r)


# This function returns the amount of a that b's projection on a covers (unit-less)
def component(a, b):
    return dotProduct(a, b) / dotProduct(a, a)


# This function returns the magnitude of the vector a
def magnitude(a):
    return math.sqrt((a[0]) ** 2 + (a[1]) ** 2 + (a[2]) ** 2)


# Determines if the lattice vectors are linearly independent
def latticeValid(a, b, c):
    M = Matrix([[a[0], a[1], a[2]], [b[0], b[1], b[2]], [c[0], c[1], c[2]]])
    return M.det() != 0


# This function returns the lattice vectors given the file exists within the same directory and the lattice vectors
# are linearly independent
def form(filename):
    # Check to ensure geometry.in file exists
    if not os.path.exists(filename):
        print("ERROR: Given geometry.in file does not exist")
        return

    # Reads in lattice vector information from geometry.in file
    with open(filename, "r") as f:
        lv = []
        for ln in f:
            if ln.startswith("lattice_vector"):
                lv.append(ln[15:])

    # Check to ensure exactly 3 lattice vectors were found
    if len(lv) != 3:
        print("ERROR: Wrong number of lattice vectors in input geometry.in file")
        return

    # Formats lattice vectors into separate variables
    aS = lv[0].split()
    bS = lv[1].split()
    cS = lv[2].split()
    a = float(aS[0]), float(aS[1]), float(aS[2])
    b = float(bS[0]), float(bS[1]), float(bS[2])
    c = float(cS[0]), float(cS[1]), float(cS[2])

    if not latticeValid(a, b, c):
        print("ERROR: Lattice vectors in geometry.in file linearly dependent")
        return
    return a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]


# Main function, input filepath of geometry.in file as well as the vector (t1,t2,t3) expressed in lattice real-space
# of the desired band direction
def generate(filename, target):
    # The form function returns a 9-element tuple containing the lattice vectors from the given geometry.in file
    lattice = form(filename)
    a = (lattice[0], lattice[1], lattice[2])
    b = (lattice[3], lattice[4], lattice[5])
    c = (lattice[6], lattice[7], lattice[8])

    # The vector (A1,A2,A3) is the target that the result will be parallel to
    A1 = (target[0] * a[0] + target[1] * b[0] + target[2] * c[0])
    A2 = (target[0] * a[1] + target[1] * b[1] + target[2] * c[1])
    A3 = (target[0] * a[2] + target[1] * b[2] + target[2] * c[2])
    target_vector = (A1, A2, A3)
    print("target: ", target_vector)

    # B1, B2, and B3 define the three reciprocal unit cell vectors. V is the volume of the unit cell.
    V = abs(dotProduct(a, crossProduct(b, c)))
    B1 = tuple((2 * math.pi / V) * x for x in crossProduct(b, c))
    B2 = tuple((2 * math.pi / V) * x for x in crossProduct(c, a))
    B3 = tuple((2 * math.pi / V) * x for x in crossProduct(a, b))

    # M defines the augmented matrix consisting of B1, B2, and B3 as columns with (A1,A2,A3) as the augmented column.
    # We know this is full rank by definition. Thus, row reducing produces the desired result
    M = Matrix([[B1[0], B2[0], B3[0], A1], [B1[1], B2[1], B3[1], A2], [B1[2], B2[2], B3[2], A3]])
    M_rref = M.rref()
    result = (M_rref[0][3], M_rref[0][7], M_rref[0][11])
    print("initial f vals: ", result)
    final = scale(target_vector, result, B1, B2, B3)
    print("f vals are: " + str(final))
    return final


f = input("Enter filepath to desired geometry.in file: ")
t1 = input("Enter x Coordinate in real lattice space: ")
t2 = input("Enter y Coordinate in real lattice space: ")
t3 = input("Enter z Coordinate in real lattice space: ")
generate(f, (float(t1), float(t2), float(t3)))