import math
import os
import numpy as np
import Rotation


# This function writes constrain relaxation flags on all types of atoms specified in the atom field using the flag
# passed in
def constrain(filepath, target, atoms, flag):
    with open(filepath, "r") as f:
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln)
            elif ln.startswith("lattice"):
                st.append(ln)
    with open(target, "w") as f:
        for x in st:
            f.write(x)
        for i in lv:
            if i.startswith("constrain") or i.startswith("initial"):
                f.write(i)
            else:
                f.write(i)
                temp = i.split()
                for j in atoms:
                    if j == temp[3]:
                        f.write("constrain_relaxation " + flag + "\n")


# This function outputs the spin splitting observed in the given band file, set to true for spin splitting
def bandProcess(filepath: object, boo, bandlength) -> object:
    min = 5.0
    prev = 0
    minVal = "idk"
    with open(filepath, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)
        j = 1
        for x in lv:
            #print("Start of k-point", j)
            temp = x.split()
            i = 0
            for y in temp:
                if 0 < i < 4:
                    p = 1
                    #print("coord: ", y)
                if (i % 2) == 0 and y == "0.00000":
                    start = i - 7
                    break
                i += 1
            i = 0
            for y in temp:
                if start < i < start + 13:
                    if (i % 2) == 0:
                        t = y
                    else:
                        p = 0
                        #print(t,y)
                    if i == start + 8:
                        if float(y) < float(min):
                            min = y
                            minVal = str(j)
                            prev = float(temp[start + 10])
                        elif (float(y) == float(min)):
                            if (float(temp[start + 10]) > prev):
                                min = y
                                minVal = str(j)
                                prev = float(temp[start + 10])
                i += 1
            j += 1
    #print("min of",str(min),"found at",minVal)
    #print("previous",prev)
    if boo:
        print(str(float(prev) - float(min)))
    else:
        if float(minVal) >= 150:
            minVal = 300 - float(minVal)
            working = float(minVal) / 300
            working *= bandlength
            print(working)
        elif float(minVal) < 150:
            working = float(minVal) / 300
            working *= bandlength
            print(working)
        else:
            print("sumn strange afoot")


def bandGap(filepath: object) -> object:
    min = 5.0
    max = 0.0
    minVal = "idk"
    maxVal = "idk either"
    with open(filepath, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)
        j = 1
        for x in lv:
            print("Start of k-point", j)
            temp = x.split()
            i = 0
            for y in temp:
                if 0 < i < 4:
                    p = 1
                    print("coord: ", y)
                if (i % 2) == 0 and y == "0.00000":
                    start = i - 7
                    break
                i += 1
            i = 0
            for y in temp:
                if start < i < start + 17:
                    if (i % 2) == 0:
                        t = y
                    else:
                        p = 0
                        print(t,y)
                    if i == start + 8:
                        if float(y) < float(min):
                            min = y
                            minVal = str(j)
                            prev = temp[start + 10]
                    if i == start + 14:
                        if float(y) > float(max):
                            max = y
                            maxVal = str(j)
                i += 1
            j += 1
    print("min of",str(min),"found at",minVal)
    print("max of",str(max),"found at",maxVal)
    print("previous",prev)
    print("diff", str(float(min)-float(prev)))
    print(float(max) - float(min))


# This function combines two geometry files which share a common atom by centering both files around that atom and
# then writing all the atoms to a new file
def combine(file1, file2, writefile, index1, index2):
    Rotation.recenter(file1, file1, index1, 0, 0, 0)
    Rotation.recenter(file2, file2, index2, 0, 0, 0)
    with open(file1, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)
    with open(file2, "r") as f:
        i = 1
        for ln in f:
            if not (i == index2):
                lv.append(ln)
            i += 1
    with open(writefile, "w") as f:
        f.writelines(lv)


# This function writes lattice vectors to a file scaled based on a distortion factor caused by in plane distortions
def addLattice(filename, writefile, original, radians):
    with open(filename, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)
    new = original * abs(math.cos(radians))
    other = 0 * abs(math.cos(radians))
    with open(writefile, "w") as f:
        f.write("lattice_vector " + str(new) + " " + str(other) + " " + "0" + "\n")
        f.write("lattice_vector " + str(other) + " " + str(new) + " 0" + "\n")
        f.write("lattice_vector 0 0 50 \n")
        f.writelines(lv)


# This function rotates a file and then adds it to a given file
def rotAdd(file1, file2, radians, axis, node1, node2, target):
    rotFile = str(file2 + "rot")
    Rotation.rotate(file2, rotFile, radians, axis)
    combine(file1, rotFile, target, node1, node2)
    # addLattice(target, target, 11.77316030, radians / 2)
    os.remove(rotFile)


def addCsInit(file, writefile):
    with open(file, "r") as f:
        lv = []
        i = -2
        for ln in f:
            lv.append(ln)
            if i == 2:
                temp = ln.split()
                x1 = temp[2]
            if i == 13:
                temp = ln.split()
                x2 = temp[2]
            if i == 6:
                temp = ln.split()
                x3 = temp[2]
            if i == 19:
                temp = ln.split()
                x4 = temp[2]
            if i == 11:
                temp = ln.split()
                y1 = temp[1]
            if i == 5:
                temp = ln.split()
                y2 = temp[1]
            if i == 8:
                temp = ln.split()
                y3 = temp[1]
            if i == 1:
                temp = ln.split()
                y4 = temp[1]
            i += 1
    with open(writefile, "w") as f:
        f.writelines(lv)
        print(x1)
        print(y1)
        x_one = (float(x1) + float(x2)) / 2
        x_two = (float(x3) + float(x4)) / 2
        y_one = (float(y1) + float(y2)) / 2
        y_two = (float(y3) + float(y4)) / 2
        f.write("atom " + str(y_one) + " " + str(x_one) + " 2.57797604083 Cs\n")
        f.write("atom " + str(y_one) + " " + str(x_two) + " 2.57797604083 Cs\n")
        f.write("atom " + str(y_two) + " " + str(x_one) + " 2.57797604083 Cs\n")
        f.write("atom " + str(y_two) + " " + str(x_two) + " 2.57797604083 Cs\n")
        f.write("atom " + str(y_one) + " " + str(x_one) + " -2.57797604417 Cs\n")
        f.write("atom " + str(y_one) + " " + str(x_two) + " -2.57797604417 Cs\n")
        f.write("atom " + str(y_two) + " " + str(x_one) + " -2.57797604417 Cs\n")
        f.write("atom " + str(y_two) + " " + str(x_two) + " -2.57797604417 Cs\n")


def addCsOOP(file, writefile, val):
    with open(file, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)

    with open(writefile, "w") as f:
        f.writelines(lv)
        f.write("atom " + str(-val) + " " + str(val) + " 2.56524642500000 Cs\n")
        f.write("atom " + str(-val) + " " + str(val) + " -2.56524642500000 Cs\n")
        f.write("atom " + str(-val) + " " + str(-val) + " 2.56524642500000 Cs\n")
        f.write("atom " + str(-val) + " " + str(-val) + " -2.56524642500000 Cs\n")
        f.write("atom " + str(val) + " " + str(val) + " 2.56524642500000 Cs\n")
        f.write("atom " + str(val) + " " + str(val) + " -2.56524642500000 Cs\n")
        f.write("atom " + str(val) + " " + str(-val) + " 2.56524642500000 Cs\n")
        f.write("atom " + str(val) + " " + str(-val) + " -2.56524642500000 Cs\n")


def betweenTotal(p1, p2, p3):
    v1 = p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]
    v2 = p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)


def between(p1, p2, p3):
    v1 = p1[0] - p2[0], p1[1] - p2[1]
    v2 = p3[0] - p2[0], p3[1] - p2[1]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    # print(u1, u2)
    # print("dp", dot_product)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)


def betweenOutV(p1, p2, p3):
    v1 = p1[1] - p2[1], p1[2] - p2[2]
    v2 = p3[1] - p2[1], p3[2] - p2[2]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)


def betweenOutH(p1, p2, p3):
    v1 = p1[0] - p2[0], p1[2] - p2[2]
    v2 = p3[0] - p2[0], p3[2] - p2[2]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)


# Alters the position of p2 to be a specific angle with p1 and p3. Varies in the direction specified by vert
def angle(filename, target, p1, p2, p3, base, angle, full, vert):
    # note there is a known bug where if calculating where angle > base but angle < 180,
    # it will find the solution where angle = 360 - angle
    # one workaround is just slightly messing with the delta initial value,
    # this often works
    with open(filename, "r") as f:
        lats = []
        at = []
        for ln in f:
            if ln.startswith("lattice"):
                lats.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(target, "w") as f:
        f.writelines(lats)
        f.write("\n")
        i = 1
        for x in at:  # This for loop pulls the atoms that the angle is being changed between
            if i == p1:
                t1 = x.split()
            elif i == p2:
                t2 = x.split()
            elif i == p3:
                t3 = x.split()
            i += 1

        if vert:
            pt1 = float(t1[2]), float(t1[1])  # is the x or y coordinate that is being used for angle calculations
            pt2 = float(t2[2]), float(t2[1])
            lat = lats[1].split()  # to be used in the case of "non-full" calculations
            if full:
                pt3 = float(t3[2]), float(t3[1])
                down = True
            else:
                pt3 = float(t3[2]) - float(lat[2]), float(t3[1]) - float(lat[1])
                down = True
        else:
            pt1 = float(t1[1]), float(t1[2])
            pt2 = float(t2[1]), float(t2[2])
            lat = lats[0].split()
            if full:
                pt3 = float(t3[1]), float(t3[2])
                down = False
            else:
                pt3 = float(t3[1]) - float(lat[1]), float(t3[2]) - float(lat[2])
                down = False

        diff = 10
        delta = 1

        limit = False
        caution = False
        if angle > base:
            caution = True
            if down:
                down = False
            else:
                down = True
            if angle > 180:
                angle = 360 - angle
                limit = True
                print("Limit = True")

        if down and (not caution or limit):
            pt2 = pt2[0], pt2[1] - delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        elif not down and (not caution or limit):
            pt2 = pt2[0], pt2[1] + delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        delta = delta / 2

        count = 1

        debug = False

        div = 2

        if debug:
            print(pt1)
            print(down)

        if limit:
            while diff > 0.0000000000005:
                delta = delta / 1.01

                if debug:
                    print()
                    print(diff)
                    print(delta)
                    print(pt2)
                count += 1

                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta

                if count == 3000:
                    print("broke due to count == 30000")
                    break

                a1 = abs(between(pt1, temp1, pt3) - angle)
                a2 = abs(between(pt1, temp2, pt3) - angle)

                if down:
                    if temp1[1] > pt1[1]:
                        continue
                    else:
                        if a1 < a2:
                            diff = a1
                            pt2 = temp1
                            if abs(pt2[1] - pt1[1]) < 0.1:
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] - delta
                                diff = abs(between(pt1, pt2, pt3) - angle)
                        else:
                            diff = a2
                            pt2 = temp2

                elif not down:
                    if temp2[1] < pt1[1]:
                        continue
                    else:
                        if a1 < a2:
                            diff = a1
                            pt2 = temp1
                        else:
                            diff = a2
                            pt2 = temp2
                            if abs(pt2[1] - pt1[1]) < 0.1:
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] + delta
                                diff = abs(between(pt1, pt2, pt3) - angle)

                else:
                    print("bruh moment")

        else:
            while diff > 0.00000000005:

                if count == 200000:
                    print("breaking from max steps reached")
                    break

                if (count % 1000) == 0:
                    delta = 0.5
                    div = 1.5

                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta

                a1 = abs(between(pt1, temp1, pt3) - angle)
                a2 = abs(between(pt1, temp2, pt3) - angle)
                if debug:
                    print()
                    print("diff:", diff)
                    print("delta:", delta)
                    print("temp1:", temp1)
                    print("temp2:", temp2)
                    print("a1:", a1)
                    print("a2:", a2)

                if a1 < a2 and a1 < diff:
                    if not caution or (down and temp1[1] >= pt1[1] - 0.1) or (not down and temp1[1] <= pt1[1] + 0.1):
                        diff = a1
                        pt2 = temp1
                elif a2 < diff:
                    if not caution or (not down and temp2[1] <= pt1[1] + 0.1) or (down and temp2[1] >= pt1[1] - 0.1):
                        diff = a2
                        pt2 = temp2
                delta = delta / div
                count += 1

        print(between(pt1, pt2, pt3))
        i = 1
        for x in at:  # This for prints out the atoms in the same order they were originally in
            if i == p2:
                if vert:
                    f.write("atom " + str(pt2[1]) + " " + str(pt2[0]) + " " + str(t2[3]) + " " + str(t2[4]) + "\n")
                else:
                    f.write("atom " + str(pt2[0]) + " " + str(pt2[1]) + " " + str(t2[3]) + " " + str(t2[4]) + "\n")
            else:
                f.write(x)
            i += 1


def angleOOP(filename, target, p1, p2, p3, base, angle, full, vert):
    with open(filename, "r") as f:
        lats = []
        at = []
        for ln in f:
            if ln.startswith("lattice"):
                lats.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(target, "w") as f:
        f.writelines(lats)
        f.write("\n")
        i = 1
        for x in at:  # This for loop pulls the atoms that the angle is being changed between
            if i == p1:
                t1 = x.split()
            elif i == p2:
                t2 = x.split()
            elif i == p3:
                t3 = x.split()
            i += 1

        if vert:
            pt1 = float(t1[2]), float(t1[3])  # is the x or y coordinate that is being used for angle calculations
            pt2 = float(t2[2]), float(t2[3])
            lat = lats[1].split()  # to be used in the case of "non-full" calculations
            if full:
                pt3 = float(t3[2]), float(t3[3])
                down = False
            else:
                pt3 = float(t3[2]) - float(lat[2]), float(t3[3]) - float(lat[3])
                down = False
        else:
            pt1 = float(t1[1]), float(t1[3])
            pt2 = float(t2[1]), float(t2[3])
            lat = lats[0].split()
            if full:
                pt3 = float(t3[1]), float(t3[3])
                down = True
            else:
                pt3 = float(t3[1]) - float(lat[1]), float(t3[3]) - float(lat[3])
                down = True

        diff = 10
        delta = 1

        limit = False
        caution = False
        if angle > base:
            caution = True
            if down:
                down = False
            else:
                down = True
            if angle > 180:
                angle = 360 - angle
                limit = True
                print("Limit = True")

        if down and (not caution or limit):
            pt2 = pt2[0], pt2[1] - delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        elif not down and (not caution or limit):
            pt2 = pt2[0], pt2[1] + delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        delta = delta / 2

        count = 1
        div = 2
        debug = False

        if debug:
            print("pt1:", pt1)
            print("pt2:", pt2)
            print("pt3:", pt3)
            print("down?", down)

        if limit:
            while diff > 0.0000000000005:
                delta = delta / 1.01

                if debug:
                    print()
                    print(diff)
                    print(delta)
                    print(pt2)
                count += 1

                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta

                if count == 30000:
                    print("broke due to count == 30000")
                    break

                a1 = abs(between(pt1, temp1, pt3) - angle)
                a2 = abs(between(pt1, temp2, pt3) - angle)

                if down:
                    if temp1[1] > pt1[1]:
                        continue
                    else:
                        if a1 < a2:
                            diff = a1
                            pt2 = temp1
                            if abs(pt2[1] - pt1[1]) < 0.1:
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] - delta
                                diff = abs(between(pt1, pt2, pt3) - angle)
                        else:
                            diff = a2
                            pt2 = temp2

                elif not down:
                    if temp2[1] < pt1[1]:
                        continue
                    else:
                        if a1 < a2:
                            diff = a1
                            pt2 = temp1
                        else:
                            diff = a2
                            pt2 = temp2
                            if abs(pt2[1] - pt1[1]) < 0.1:
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] + delta
                                diff = abs(between(pt1, pt2, pt3) - angle)

                else:
                    print("bruh moment")

        else:
            while diff > 0.00000000005:

                if count == 20000:
                    print("breaking from max steps reached")
                    break

                if (count % 1000) == 0:
                    delta = 0.5
                    div = 1.5

                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta

                a1 = abs(between(pt1, temp1, pt3) - angle)
                a2 = abs(between(pt1, temp2, pt3) - angle)
                if debug:
                    print()
                    print("diff:", diff)
                    print("delta:", delta)
                    print("temp1:", temp1)
                    print("temp2:", temp2)
                    print("a1:", a1)
                    print("a2:", a2)

                if a1 < a2 and a1 < diff:
                    if not caution or (down and temp1[1] >= pt1[1] - 0.1) or (not down and temp1[1] <= pt1[1] + 0.1):
                        diff = a1
                        pt2 = temp1
                elif a2 < diff:
                    if not caution or (not down and temp2[1] <= pt1[1] + 0.1) or (down and temp2[1] >= pt1[1] - 0.1):
                        diff = a2
                        pt2 = temp2
                delta = delta / div
                count += 1

        print(between(pt1, pt2, pt3))
        i = 1
        for x in at:  # This for prints out the atoms in the same order they were originally in
            if i == p2:
                if vert:
                    f.write("atom " + str(t2[1]) + " " + str(pt2[0]) + " " + str(pt2[1]) + " " + str(t2[4]) + "\n")
                else:
                    f.write("atom " + str(pt2[0]) + " " + str(t2[2]) + " " + str(pt2[1]) + " " + str(t2[4]) + "\n")
            else:
                f.write(x)
            i += 1


def charge(filename, target, atoms, charges):
    with open(filename, "r") as f:
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
            elif ln.startswith("constrain") or ln.startswith("initial") or ln.startswith("lattice"):
                lv.append(ln)
    with open(target, "w") as f:
        for i in lv:
            if i.startswith("constrain") or i.startswith("initial") or i.startswith("lattice"):
                f.write(i)
            else:
                temp = i.split()
                p1 = float(temp[0]), float(temp[1]), float(temp[2])
                final = p1
                f.write("atom " + str(final[0]) + " " + str(final[1]) + " " + str(final[2]) + " " + str(temp[3]) + "\n")
                x = 0
                for j in atoms:
                    if j == temp[3]:
                        f.write("initial_charge " + str(charges[x]) + "\n")
                    x += 1


# outputs delta beta, delta beta in plane, and delta beta out of plane in that order
def angleInfo(filename, p1, p2, p3, vert):
    with open(filename, "r") as f:
        lv = [0, 0, 0]
        i = 1
        for ln in f:
            if ln.startswith("atom"):
                if (i == p1):
                    lv[0] = ln
                if (i == p2):
                    lv[1] = ln
                if (i == p3):
                    lv[2] = ln
                i += 1
        t1 = lv[0].split()
        t2 = lv[1].split()
        t3 = lv[2].split()
        a1 = float(t1[1]), float(t1[2]), float(t1[3])
        a2 = float(t2[1]), float(t2[2]), float(t2[3])
        a3 = float(t3[1]), float(t3[2]), float(t3[3])
        print(betweenTotal(a1, a2, a3))
        print(between(a1, a2, a3))
        if vert:
            print(betweenOutV(a1, a2, a3))
        else:
            print(betweenOutH(a1, a2, a3))


def printBands(filestart, bandnum, bandlength):
    a = 180
    for i in range(1):
        print(a)
        print("band gap")
        for i in range(9):
            fileplus = filestart + str(a) + "/" + str(a) + "_" + str((a - 20) + (5 * i))
            fileplus += "/band100" + str(bandnum) + ".out"
            bandProcess(fileplus, True, bandlength)
        print()  # k0 value
        for i in range(9):
            fileplus = filestart + str(a) + "/" + str(a) + "_" + str((a - 20) + (5 * i))
            fileplus += "/band100" + str(bandnum) + ".out"
            bandProcess(fileplus, False, bandlength)
        print()  # band width value
        for i in range(9):
            fileplus = filestart + str(a) + "/" + str(a) + "_" + str((a - 20) + (5 * i))
            fileplus += "/band100" + str(bandnum) + ".out"
            bandGap(fileplus)
        a -= 10


def make_I_vert(readfile, writefile, M_height, Bi_height, Cs_height):
    with open(readfile, "r") as f:
        lv = []
        at = []
        for ln in f:
            if ln.startswith("lattice"):
                lv.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(writefile, "w") as f:
        f.writelines(lv)
        # temp = lv[0].split()
        # Cs_loc = float(temp[1]) / 4
        for i in at:
            temp = i.split()
            if (2.5 > float(temp[3]) > -2.5) or temp[4] == "Cs":
                f.write(i)
            if temp[4] == "Tl":
                temp[4] = "I"
                temp[3] = str(float(temp[3]) + M_height)
                f.write(' '.join(temp))
                f.write("\n")
                temp[3] = str(float(temp[3]) - 2 * M_height)
                f.write(' '.join(temp))
                f.write("\n")
            if temp[4] == "Bi":
                temp[4] = "I"
                temp[3] = str(float(temp[3]) + Bi_height)
                f.write(' '.join(temp))
                f.write("\n")
                temp[3] = str(float(temp[3]) - 2 * Bi_height)
                f.write(' '.join(temp))
                f.write("\n")
        # Cs = ["atom", str(Cs_loc), str(Cs_loc), str(Cs_height), "Cs"]
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[3] = str(-float(Cs[3]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[1] = str(-float(Cs[1]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[3] = str(-float(Cs[3]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[1] = str(-float(Cs[1]))
        # Cs[2] = str(-float(Cs[2]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[3] = str(-float(Cs[3]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[1] = str(-float(Cs[1]))
        # f.write(' '.join(Cs))
        # f.write("\n")
        # Cs[3] = str(-float(Cs[3]))
        # f.write(' '.join(Cs))
        # f.write("\n")


def addCs(readfile, writefile, Cs_height):
    with open(readfile, "r") as f:
        lv = []
        at = []
        for ln in f:
            if ln.startswith("lattice"):
                lv.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(writefile, "w") as f:
        f.writelines(lv)
        for i in at:
            temp = i.split()
            if temp[4] != "Cs":
                f.write(i)
                f.write("\n")
        temp = lv[0].split()
        Cs_loc = float(temp[1]) / 4
        Cs = ["atom", str(Cs_loc), str(Cs_loc), str(Cs_height), "Cs"]
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[3] = str(-float(Cs[3]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[1] = str(-float(Cs[1]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[3] = str(-float(Cs[3]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[1] = str(-float(Cs[1]))
        Cs[2] = str(-float(Cs[2]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[3] = str(-float(Cs[3]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[1] = str(-float(Cs[1]))
        f.write(' '.join(Cs))
        f.write("\n")
        Cs[3] = str(-float(Cs[3]))
        f.write(' '.join(Cs))
        f.write("\n")


def delta_d_calc(readfile, center, x1, x2, y1, y2, z1, z2, full):
    with open(readfile, "r") as f:
        lv = []
        at = []
        extra = 0
        i = 1
        for ln in f:
            if ln.startswith("lattice_vector"):
                lv.append(ln)
            if ln.startswith("atom"):
                if i == center:
                    cent = ln
                elif (i == x1) or (i == x2) or (i == y1) or (i == z1) or (i == z2):
                    at.append(ln)
                elif i == y2:
                    extra = ln
                i += 1
        if not full:
            lat = lv[1].split()
            extra = extra.split()
            extra[1] = str(float(extra[1]) + float(lat[1]))
            extra[2] = str(float(extra[2]) + float(lat[2]))
            extra[3] = str(float(extra[3]) + float(lat[3]))
            extra = ' '.join(extra)
        at.append(extra)

        distances = []
        # print(at)
        for x in at:
            distances.append(distance(cent, x))

        avg = 0
        for x in distances:
            avg += x
        avg = avg / 6
        avg2 = avg ** 2

        result = 0
        for x in distances:
            result += ((x - avg) ** 2) / avg2

        print(result / 6)


def old_sigma_square_calc(readfile, center, x1, x2, y1, y2, z1, z2, full):
    with open(readfile, "r") as f:
        lv = []
        at = []
        extra = 0
        i = 1
        for ln in f:
            if ln.startswith("lattice_vector"):
                lv.append(ln)
            if ln.startswith("atom"):
                if i == center:
                    center = ln
                elif i == x1:
                    x1 = ln
                elif i == x2:
                    x2 = ln
                elif i == y1:
                    y1 = ln
                elif i == y2:
                    y2 = ln
                elif i == z1:
                    z1 = ln
                elif i == z2:
                    z2 = ln
                i += 1
        if not full:
            lat = lv[1].split()
            y2 = y2.split()
            y2[1] = str(float(y2[1]) + float(lat[1]))
            y2[2] = str(float(y2[2]) + float(lat[2]))
            y2[3] = str(float(y2[3]) + float(lat[3]))
            y2 = ' '.join(y2)

        angles = []
        angles.append(tri_angle(x1, center, y1))
        angles.append(tri_angle(x1, center, y2))
        angles.append(tri_angle(x1, center, z1))
        angles.append(tri_angle(x1, center, z2))

        angles.append(tri_angle(y2, center, x2))
        angles.append(tri_angle(y2, center, z1))
        angles.append(tri_angle(y2, center, z2))

        angles.append(tri_angle(x2, center, y1))
        angles.append(tri_angle(x2, center, z1))
        angles.append(tri_angle(x2, center, z2))

        angles.append(tri_angle(y1, center, z1))
        angles.append(tri_angle(y1, center, z2))

        avg = 0
        for x in angles:
            avg += x
        avg = avg / 12

        result = 0
        for x in angles:
            result += ((x - avg) ** 2)

        print(result / 12)


def sigma_square_calc(readfile, xmin, xmax, ymin, ymax, center, bismuth):
    with open(readfile, "r") as f:
        at = []
        lv = []
        for ln in f:
            if ln.startswith("atom"):
                at.append(ln)
            elif ln.startswith("lattice"):
                lv.append(ln)

    oct = []
    for a in at:
        temp = a.split()
        if (xmin <= float(temp[1]) <= xmax) and (ymin <= float(temp[2]) <= ymax) and temp[4] == "I":
            oct.append(a)
        elif (xmin <= float(temp[1]) <= xmax) and (ymin <= float(temp[2]) <= ymax) and temp[4] == center:
            center = a
        if (-1 <= float(temp[1]) <= 1) and (float(temp[2]) <= -1):
            possible = a.split()
    if bismuth:
        tempLV = lv[1].split()
        new = (float(possible[1]) + float(tempLV[1]), float(possible[2]) + float(tempLV[2]),
               float(possible[3]) + float(tempLV[3]))
        temp = "atom " + str(new[0]) + " " + str(new[1]) + " " + str(new[2]) + " I"
        oct.append(temp)

    angles = []
    i = 0
    for x in oct:
        j = 0
        for y in oct:
            if i < j:
                temp = tri_angle(x, center, y)
                if 45 < temp < 135:
                    angles.append(temp)
            j += 1
        i += 1

    if len(angles) != 12:
        print("SAD!", len(angles))
        print(angles)

    result = 0
    for x in angles:
        result += ((x - 90) ** 2)

    print(result / 12)


def tri_angle(t1, t2, t3):
    t1 = t1.split()
    t2 = t2.split()
    t3 = t3.split()

    p1 = float(t1[1]), float(t1[2]), float(t1[3])
    p2 = float(t2[1]), float(t2[2]), float(t2[3])
    p3 = float(t3[1]), float(t3[2]), float(t3[3])

    v1 = p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]
    v2 = p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    # print(u1, u2)
    # print("dp", dot_product)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)


def distance(a, b):
    a = a.split()
    b = b.split()
    a_num = float(a[1]), float(a[2]), float(a[3])
    b_num = float(b[1]), float(b[2]), float(b[3])
    result = ((a_num[0] - b_num[0]) ** 2) + ((a_num[1] - b_num[1]) ** 2) + ((a_num[2] - b_num[2]) ** 2)
    return math.sqrt(result)


def centroid(readfile, xmin, xmax, ymin, ymax, center, bismuth):
    with open(readfile, "r") as f:
        at = []
        lv = []
        for ln in f:
            if ln.startswith("atom"):
                at.append(ln)
            elif ln.startswith("lattice"):
                lv.append(ln)
    oct = []
    for a in at:
        temp = a.split()
        if (xmin <= float(temp[1]) <= xmax) and (ymin <= float(temp[2]) <= ymax) and temp[4] == "I":
            oct.append(a)
        elif (xmin <= float(temp[1]) <= xmax) and (ymin <= float(temp[2]) <= ymax) and temp[4] == center:
            center = a.split()
        if (-1 <= float(temp[1]) <= 1) and (float(temp[2]) <= -1):
            possible = a.split()
    if bismuth:
        tempLV = lv[1].split()
        new = (float(possible[1]) + float(tempLV[1]), float(possible[2]) + float(tempLV[2]),
               float(possible[3]) + float(tempLV[3]))
        temp = "atom " + str(new[0]) + " " + str(new[1]) + " " + str(new[2]) + " I"
        oct.append(temp)
    x = 0.0
    y = 0.0
    z = 0.0
    for o in oct:
        # print(o)
        temp = o.split()
        x += float(temp[1])
        y += float(temp[2])
        z += float(temp[3])
    # print(center)
    x = float(center[1]) - (x / 6)
    y = float(center[2]) - (y / 6)
    z = float(center[3]) - (z / 6)
    output = str(x) + ", " + str(y) + ", " + str(z)
    if len(oct) != 6:
        print("AHHHHHHHHHHHHHHHHH", len(oct))
    print(output)


def flipYandZ(readfile, writefile):
    with open(readfile, "r") as f:
        at = []
        lv = []
        for ln in f:
            if ln.startswith("atom"):
                at.append(ln)
            elif ln.startswith("lattice"):
                lv.append(ln)
    with open(writefile, "w") as f:
        for l in lv:
            temp = l.split()
            outp = "lattice_vector\t" + temp[1] + "\t" + temp[3] + "\t" + temp[2] + "\n"
            f.write(outp)
        for a in at:
            temp = a.split()
            outp = "atom\t" + temp[1] + "\t" + temp[3] + "\t" + temp[2] + "\t" + temp[4] + "\n"
            f.write(outp)
        print("Flipped Y and Z :)")