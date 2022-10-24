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


# This function outputs the spin splitting observed in the given band file
def bandProcess(filepath: object, boo) -> object:
    min = 5.0
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
                    start = i-7
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
                            prev = temp[start+10]
                i += 1
            j += 1
    #print("min of",str(min),"found at",minVal)
    #print("previous",prev)
    if boo:
        print(str(float(prev)-float(min)))
    else:
        print(minVal)

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
            #print("Start of k-point", j)
            temp = x.split()
            i = 0
            for y in temp:
                if 0 < i < 4:
                    p = 1
                    #print("coord: ", y)
                if (i % 2) == 0 and y == "0.00000":
                    start = i-7
                    break
                i += 1
            i = 0
            for y in temp:
                if start < i < start + 17:
                    if (i % 2) == 0:
                        t = y
                    else:
                        p = 0
                        #print(t,y)
                    if i == start + 8:
                        if float(y) < float(min):
                            min = y
                            minVal = str(j)
                            prev = temp[start+10]
                    if i == start + 14:
                        if float(y) > float(max):
                            max = y
                            maxVal = str(j)
                i += 1
            j += 1
    #print("min of",str(min),"found at",minVal)
    #print("max of",str(max),"found at",maxVal)
    #print("previous",prev)
    #print("diff", str(float(min)-float(prev)))
    print(float(max)-float(min))

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
        f.write("lattice_vector "+str(new)+" "+str(other)+" "+"0"+"\n")
        f.write("lattice_vector "+str(other)+" "+str(new)+" 0"+"\n")
        f.write("lattice_vector 0 0 50 \n")
        f.writelines(lv)

# This function rotates a file and then adds it to a given file
def rotAdd(file1, file2, radians, axis, node1, node2, target):
    rotFile = str(file2+"rot")
    Rotation.rotate(file2, rotFile, radians, axis)
    combine(file1, rotFile, target, node1, node2)
    #addLattice(target, target, 11.77316030, radians / 2)
    os.remove(rotFile)


def addCs(file,writefile):
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
        x_one = (float(x1) + float(x2))/2
        x_two = (float(x3) + float(x4))/2
        y_one = (float(y1) + float(y2))/2
        y_two = (float(y3) + float(y4))/2
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
    return (np.arccos(dot_product))*360/(2*math.pi)

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

# Alters the position of p2 to be a specific angle with p1 and p3. Varies in the coord direction
def angle(filename, target, p1, p2, p3, full, angle, coord):
    with open(filename, "r") as f:
        lv = []
        for ln in f:
            lv.append(ln)
    with open(target, "w") as f:
        i = -2
        l = 0
        for x in lv:
            if not(x.startswith("lattice") or x.startswith("atom")):
                continue
            if i == p1:
                t1 = x.split()
                f.write(x)
            elif x.startswith("lattice"):
                l += 1
                f.write(x)
                if coord == 1 and l == 1:
                    lat = x.split()
                if coord == 2 and l == 2:
                    lat = x.split()
                if coord == 3 and l == 3:
                    lat = x.split()
            elif i == p2:
                t2 = x.split()
            elif i == p3:
                t3 = x.split()
                f.write(x)
            else:
                f.write(x)
            i += 1
        pt1 = float(t1[1]), float(t1[2]), float(t1[3])
        pt2 = float(t2[1]), float(t2[2]), float(t2[3])
        if full:
            pt3 = float(t3[1]), float(t3[2]), float(t3[3])
        else:
            pt3 = float(t3[1]) - float(lat[1]), float(t3[2]) - float(lat[2]), float(t3[3]) - float(lat[3])
        diff = 1
        delta = .3
        y = 0
        if coord == 2 and angle < between(pt1, pt2, pt3):
            down = True
        elif coord == 1 and angle > between(pt1, pt2, pt3):
            down = True
        elif coord == 3:
            print("unfinished code, results bad")
        else:
            down = False
        limit = False
        if angle > 180:
            d = angle - 180
            angle = 180 - d
            limit = True
        while diff > 0.00000000005:
            if coord == 2:
                temp1 = pt2[0] + delta, pt2[1]
                temp2 = pt2[0] - delta, pt2[1]
            elif coord == 1:
                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta
            a1 = abs(between(pt1, temp1, pt3) - angle)
            a2 = abs(between(pt1, temp2, pt3) - angle)
            if y == 0:
                if down:
                    diff = a2
                    pt2 = temp2
                else:
                    diff = a1
                    pt2 = temp1
            else:
                if limit:
                    if not down and coord == 2:
                        if temp2[0] < pt1[0]:
                            pt2 = temp1
                            diff = a1
                        else:
                            diff = min(a1, a2)
                            if a1 < a2:
                                pt2 = temp1
                            else:
                                pt2 = temp2
                    elif down and coord == 1:
                        if temp1[1] > pt1[1]:
                            pt2 = temp2
                            diff = a2
                        else:
                            diff = min(a1, a2)
                            if a1 < a2:
                                pt2 = temp1
                            else:
                                pt2 = temp2
                    else:
                        print("bruh moment")
                else:
                    diff = min(a1, a2)
                    if a1 < a2:
                        pt2 = temp1
                    else:
                        pt2 = temp2
            y += 1
            delta = delta / 2
        print(between(pt1, pt2, pt3))
        f.write("atom " + str(pt2[0]) + " " + str(pt2[1]) + " " + str(t2[3]) + " " + str(t2[4]) + "\n")


def angleOOP(filename, target, p1, p2, p3, base, angle, full, vert):
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
        for x in at:            # This for loop pulls the atoms that the angle is being changed between
            if i == p1:         # And writes all atoms except the one that needs to move to the output file
                t1 = x.split()
                f.write(x)
            elif i == p2:
                t2 = x.split()
            elif i == p3:
                t3 = x.split()
                f.write(x)
            else:
                f.write(x)
            i += 1


        if vert:                                # This if-else clause sets up the atoms depending on whether it
            pt1 = float(t1[2]), float(t1[3])    # is the x or y coordinate that is being used for angle calculations
            pt2 = float(t2[2]), float(t2[3])    # notice it pulling the corresponding lattice vector as well
            lat = lats[1].split()               # to be used in the case of "non-full" calculations
        else:
            pt1 = float(t1[1]), float(t1[3])
            pt2 = float(t2[1]), float(t2[3])
            lat = lats[0].split()

        if full and vert:
            pt3 = float(t3[2]), float(t3[3])
            down = False
        elif full:
            pt3 = float(t3[1]), float(t3[3])
            down = True
        elif vert:
            pt3 = float(t3[2]) - float(lat[2]), float(t3[3]) - float(lat[3])
            down = False
        elif not vert:
            pt3 = float(t3[1]) - float(lat[1]), float(t3[3]) - float(lat[3])
            down = True
        else:
            print("bruh")

        diff = 1
        delta = .5

        # Does this have a use?
        limit = False
        if angle > 180:
            angle = 360 - angle
            limit = True
            print("Limit = True")
            if down:
                down = False
            else:
                down = True
        elif angle > base:
            if down:
                down = False
            else:
                down = True

        if down:
            pt2 = pt2[0], pt2[1] - delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        else:
            pt2 = pt2[0], pt2[1] + delta
            diff = abs(between(pt1, pt2, pt3) - angle)
        delta = delta / 2

        count = 0
        if limit:
            while diff > 0.00000000005:
                delta = delta / 1.01
                #print(diff)
                count += 1

                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta


                if count == 30000:
                    print("broke due to count == 3000")
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
                            if(pt2[1] > -0.1):
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] - delta
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
                            if (pt2[1] < 0.1):
                                delta = delta / 1.01
                                pt2 = pt2[0], pt2[1] + delta
                else:
                    print("bruh moment")

        else:
            while diff > 0.00000000005:
                temp1 = pt2[0], pt2[1] + delta
                temp2 = pt2[0], pt2[1] - delta

                a1 = abs(between(pt1, temp1, pt3) - angle)
                a2 = abs(between(pt1, temp2, pt3) - angle)

                if a1 < a2:
                    diff = a1
                    pt2 = temp1
                else:
                    diff = a2
                    pt2 = temp2
                delta = delta / 2

        print(between(pt1, pt2, pt3))
        if vert:
            f.write("atom " + str(t2[1]) + " " + str(pt2[0]) + " " + str(pt2[1]) + " " + str(t2[4]) + "\n")
        else:
            f.write("atom " + str(pt2[0]) + " " + str(t2[2]) + " " + str(pt2[1]) + " " + str(t2[4]) + "\n")


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
        lv = [0,0,0]
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
