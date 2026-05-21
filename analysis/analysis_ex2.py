#同時測定の角度分布(円柱)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import *

def Om(l, h, R):
    eta2 = 4 * R * h/((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt( (h + R)**2 + l**2 )

N = 100
th = np.linspace(0, 75, N)
thr = np.radians(th)

L = 5
Rs = 1.27
Rb = 2.54
c = 3
d = 0.1
r = 0.6

x = np.linspace(- d/2, d/2, N)
rho = np.linspace(0, r, N)
phi = np.linspace(0, 2 * np.pi, N)
X, RHO, PHI = np.meshgrid(x, rho, phi, indexing="ij")


#l = L + x[i] * np.cos(thr) + rho[j] * np.sin(thr) * np.sin(phi[k])
#h = np.abs(x[i] * np.sin(thr) - rho[j] * np.sin(phi[k]) * np.cos(phi[k]))
#xhp = x[i] + L * np.cos(thr)
#yhp = rho[j] * np.sin(phi[k]) + L * np.sin(thr)
#zhp = rho[j] * np.cos(phi[k])
#HP = np.sqrt(xhp**2 + yhp**2 + zhp**2)
#yT = yhp * (c + L * np.cos(thr)) / xhp - L * np.sin(thr)
#zT = zhp * (c + L * np.cos(thr)) / xhp