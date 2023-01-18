import math
import numpy as np


def vecSub(p1, p2):
    return p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]


def vecCross(p1, p2):
    return p1[1] * p2[2] - p1[2] * p2[1], p1[2] * p2[0] - p1[0] * p2[2], p1[0] * p2[1] - p1[1] * p2[0]


def vecNormalize(p1):
    mag = math.sqrt((p1[0] ** 2) + (p1[1] ** 2) + (p1[2] ** 2))
    if mag == 0:
        print("ERROR: vector magnitude 0")
        return 0, 0, 0
    return p1[0] / mag, p1[1] / mag, p1[2] / mag


def vecDot(p1, p2):
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]


def projPTonPlane(p1, n):
    n = vecNormalize(n)
    mag = vecDot(p1, n)
    diff = (n[0] * mag, n[1] * mag, n[2] * mag)
    return vecSub(p1, diff)


def pointAngle(p1, p2, p3):
    v1 = p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]
    v2 = p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(u1, u2)
    return (np.arccos(dot_product)) * 360 / (2 * math.pi)