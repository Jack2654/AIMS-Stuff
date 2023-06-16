import VectorToolkit as vt
import aimstools


def angleInfo(p1, p2, p3):
    debug = 0
    p3 = vt.vecSub(p3, p1)
    p2 = vt.vecSub(p2, p1)
    p1 = (0, 0, 0)
    v1 = vt.vecSub(p3, p1)
    t1 = (v1[0], v1[1], 0)
    v2 = vt.vecCross(t1, (0, 0, 1))
    n = vt.vecNormalize(vt.vecCross(v1, v2))
    p = vt.vecNormalize(vt.vecCross(v1, n))
    IP = vt.projPTonPlane(p2, n)
    OOP = vt.projPTonPlane(p2, p)
    a = vt.pointAngle(p1, p2, p3)
    IPa = vt.pointAngle(p1, IP, p3)
    OOPa = vt.pointAngle(p1, OOP, p3)

    if debug:
        print("p1:\t", p1)
        print("p2:\t", p2)
        print("p3:\t", p3)
        print("t1:\t", t1)
        print("v1:\t", v1)
        print("v2:\t", v2)
        print("n:\t", n)
        print("plane equation is " + str(n[0]) + "x + " + str(n[1]) + "y + " + str(n[2]) + "z = 0")
        print("plane equation is " + str(p[0]) + "x + " + str(p[1]) + "y + " + str(p[2]) + "z = 0")
        print("p2 IP:\t\t", IP)
        print("p2 OOP:\t\t", OOP)
    return a, IPa, OOPa


def manyLocAngles(source, geometry, lat, mod, plus):  # this method finds the angle info in the format angle, IP, OOP
    if (mod != "a") and (mod != "f") and (mod != "n") and (mod != "s"):
        print("Passed in bad fourth argument [mod], see commented code for correct usage")
    with open(geometry, "r") as f:  # for many angles with atoms specified in the source file
        lats = []  # the source file has one set of three points per line with each atom
        at = []  # number being whitespace delimited
        for ln in f:  # in particular angle calculations use a locally defined
            if ln.startswith("lattice"):  # of IP and OOP
                lats.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(source, "r") as f:
        pts = []
        for ln in f:
            pts.append(ln)
    lat0 = lats[0].split()
    lat1 = lats[1].split()
    lat2 = lats[2].split()
    lvals = float(lat0[1]), float(lat0[2]), float(lat0[3]), float(lat1[1]), float(lat1[2]), float(lat1[3]), float(
        lat2[1]), float(lat2[2]), float(lat2[3])
    i = 0
    for x in pts:
        temp = x.split()
        t1 = at[int(temp[0]) - 1].strip().split()
        t2 = at[int(temp[1]) - 1].strip().split()
        t3 = at[int(temp[2]) - 1].strip().split()
        p1 = [float(t1[1]), float(t1[2]), float(t1[3])]
        p2 = [float(t2[1]), float(t2[2]), float(t2[3])]
        p3 = [float(t3[1]), float(t3[2]), float(t3[3])]
        if (mod == "a") or (mod == "f" and (i % 2) == 0) or (mod == "s" and (i % 2) == 1):
            if plus:
                p3[0] += lvals[3 * lat]
                p3[1] += lvals[3 * lat + 1]
                p3[2] += lvals[3 * lat + 2]
            else:
                p3[0] -= lvals[3 * lat]
                p3[1] -= lvals[3 * lat + 1]
                p3[2] -= lvals[3 * lat + 2]
        a = angleInfo(p1, p2, p3)
        print(a[0], a[1], a[2])
        i += 1


def manyAbsAngles(source, geometry, lat, mod, plus):  # this method is similar to the one above except it defines IP and OOP
    with open(geometry, "r") as f:  # based on the "in-plane" lattice vectors
        lats = []  # this makes little sense in my case :(
        at = []  # mod takes the values "a" meaning all pairs need a lattice added,
        for ln in f:  # "n" meaning no pairs, "f" means every other starting with the first
            if ln.startswith("lattice"):  # and "s" means every other starting with the second
                lats.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(source, "r") as f:
        pts = []
        for ln in f:
            pts.append(ln)
    lat0 = lats[0].split()
    lat1 = lats[1].split()
    lat2 = lats[2].split()
    lvals = float(lat0[1]), float(lat0[2]), float(lat0[3]), float(lat1[1]), float(lat1[2]), float(lat1[3]), float(
        lat2[1]), float(lat2[2]), float(lat2[3])
    i = 0
    for x in pts:
        temp = x.split()
        t1 = at[int(temp[0]) - 1].strip().split()
        t2 = at[int(temp[1]) - 1].strip().split()
        t3 = at[int(temp[2]) - 1].strip().split()
        p1 = [float(t1[1]), float(t1[2]), float(t1[3])]
        p2 = [float(t2[1]), float(t2[2]), float(t2[3])]
        p3 = [float(t3[1]), float(t3[2]), float(t3[3])]
        if (mod == "a") or (mod == "f" and (i % 2) == 0) or (mod == "s" and (i % 2) == 1):
            if plus:
                p3[0] += lvals[3 * lat]
                p3[1] += lvals[3 * lat + 1]
                p3[2] += lvals[3 * lat + 2]
            else:
                p3[0] -= lvals[3 * lat]
                p3[1] -= lvals[3 * lat + 1]
                p3[2] -= lvals[3 * lat + 2]

        ang = vt.pointAngle(p1, p2, p3)
        p3[2] = 0  # this is only because I know in my case the lattice vectors make up the xy plane
        p2[2] -= p1[2]  # otherwise would actually have to do something different
        p1[2] -= p1[2]

        a = angleInfo(p1, p2, p3)
        print(ang, a[1], a[2])
        i += 1


def moveManyPoints(source, geo, lat, writefile, plus):
    with open(geo, "r") as f:
        lats = []
        at = []
        for ln in f:
            if ln.startswith("lattice"):
                lats.append(ln)
            elif ln.startswith("atom"):
                at.append(ln)
    with open(source, "r") as f:
        pts = []
        for ln in f:
            pts.append(ln.strip())
    lat0 = lats[0].split()
    lat1 = lats[1].split()
    lat2 = lats[2].split()
    lvals = float(lat0[1]), float(lat0[2]), float(lat0[3]), float(lat1[1]), float(lat1[2]), float(lat1[3]), float(
        lat2[1]), float(lat2[2]), float(lat2[3])

    with open(writefile, "w") as f:
        f.writelines(lats)
        i = 0
        for x in at:
            flag = False
            for y in pts:
                if float(y) == (i + 1):
                    flag = True
            if not flag:
                f.write(x)
            else:
                print("lalalala")
                temp = x.split()
                tempPT = [float(temp[1]), float(temp[2]), float(temp[3])]
                if plus:
                    tempPT[0] += lvals[3 * lat]
                    tempPT[1] += lvals[3 * lat + 1]
                    tempPT[2] += lvals[3 * lat + 2]
                else:
                    tempPT[0] -= lvals[3 * lat]
                    tempPT[1] -= lvals[3 * lat + 1]
                    tempPT[2] -= lvals[3 * lat + 2]
                f.write(
                    temp[0] + " " + str(tempPT[0]) + " " + str(tempPT[1]) + " " + str(tempPT[2]) + " " + temp[4] + "\n")
            i += 1


def mullikenPlotter(filepath):
    # these two lines are only necessary to make the jupyter notebooks run on binder
    # We load the respective module
    from aimstools.bandstructures import MullikenBandStructure as MBS
    # import sys
    # sys.path.insert(0, filepath)
    # We initialize the MullikenBandStructure object from a directory.
    # The bandmlkxxxx.out files can become very large, so reading them in can take up some time.
    bs = MBS(filepath, soc=True)
    import matplotlib.pyplot as plt

    # We set up a figure and some axes
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # On the first one, we draw just the band structure, in this case with spin-orbit coupling:
    contributions = ["Pb", "I"]
    colors = ["red", "blue"]
    ax1 = bs.plot(axes=ax1, color="royalblue", show_bandgap_vertices=False)
    ax1.set_title("Band Structure with SOC")
    plt.ylim([-3,3])

    # On the second one, we draw the contribution of all species overlaid:
    ax2 = bs.plot_contributions(contributions=contributions, labels=["Pb", "I"], colors=colors, mode="scatter", axes=ax2)
    ax2.set_ylabel("")
    ax2.set_yticks([])
    ax2.set_title("Projection on scatter difference")
    plt.ylim([-3, 3])

    # On the third one, we draw the difference of the contributions as a gradient:
    ax3 = bs.plot_contributions(contributions=contributions, labels=["Pb", "I"], colors=colors, mode="scatter", bands_alpha=0.4, axes=ax3)
    ax3.set_ylabel("")
    ax3.set_yticks([])
    ax3.set_title("Projection on scatter difference")
    plt.ylim([-3,3])
    plt.show()
