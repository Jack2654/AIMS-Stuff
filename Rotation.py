from scipy.spatial.transform import Rotation as R
import os.path


def rot(p1, theta, axis):
    rotation_radians = theta
    rotation_axis = axis
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    rotated_vec = rotation.apply(p1)
    return rotated_vec


# This function rotates a geometry (about zero) about the given axis using the given theta
def rotate(filepath, writepath, theta, axis):
    with open(filepath, "r") as f:
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
            elif ln.startswith("lattice"):
                st.append(ln)
    with open(writepath, "w") as f:
        f.writelines(st)
        for i in lv:
            temp = i.split()
            p1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = rot(p1, theta, axis)
            f.write("atom "+str(final[0])+" "+str(final[1])+" "+str(final[2])+" "+str(temp[3])+"\n")


def reflect(filepath):
    with open(filepath, "r") as f:
        lv = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
        for i in lv:
            temp = i.split()
            p1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = p1[0], -p1[1], p1[2]
            print("atom", final[0], final[1], final[2], temp[3])


# This function recenters the file about the given x,y,x coordinates being written to the index provided
def recenter(filepath, writepath, index, x, y, z):
    with open(filepath, "r") as f:
        j = 0
        lv = []
        st = []
        for ln in f:
            if ln.startswith("atom"):
                lv.append(ln[5:])
            else:
                st.append(ln)
        for i in lv:
            j += 1
            if j == index:
                temp = i.split()
                t1 = float(temp[0])-x, float(temp[1])-y, float(temp[2])-z
                p1 = t1
    with open(writepath, "w") as f:
        f.writelines(st)
        for i in lv:
            temp = i.split()
            t1 = float(temp[0]), float(temp[1]), float(temp[2])
            final = t1[0]-p1[0],t1[1]-p1[1],t1[2]-p1[2]
            f.write("atom " + str(final[0]) + " " + str(final[1]) + " " + str(final[2]) + " " + str(temp[3]) + "\n")