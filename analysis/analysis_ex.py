#同時測定の角度分布(点近似)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def Om(l,Rs):
    return 2 * np.pi * (1 - l / np.sqrt(Rs**2 + l**2))

th = np.linspace(0, 90, 300)
thr = th * np.pi / 360

l = 5
Rs = 1.27
Rb = 2.54
Om = Om(l, Rs)
b = np.sqrt(Om / 2 * np.pi) / np.cos(thr)
a = b / np.cos(thr)





