#同時測定の角度分布(円柱)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import *

def Om(l, h, R):
    eta2 = 4 * R * h/((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt( (h + R)**2 + l**2 )

th = np.linspace(0, 75, 300)
thr = np.radians(th)

l = 5
Rs = 1.27
Rb = 2.54
c = 3
d = 0.1
r = 0.6

x = np.linspace(- d/2, d/2, 100)
rho = np.linspace(0, r, 100)
phi = np.linspace(0, 2 * np.pi, 100)
