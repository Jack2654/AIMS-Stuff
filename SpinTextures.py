import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_state(filename):
    with open(filename, "r") as f:
        st = []
        for ln in f:
            if not ln.startswith("#"):
                st.append(ln)

    tr = []
    for x in range(len(st)):
        if x == 0:
            tr.append(st[x])
        else:
            if(st[x-1].split()[0] != st[x].split()[0]):
                tr.append(st[x])
    return tr