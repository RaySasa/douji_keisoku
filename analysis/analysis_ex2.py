#同時測定の角度分布(円柱)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import *

def Om(l, h, R):
    eta2 = 4 * R * h/((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt( (h + R)**2 + l**2 )

N = 100
thr = np.radians(np.linspace(0, 75, N))

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

# 結果保存
result = []

for th in thr:

    # 色々定義
    l = L + X * np.cos(th) + RHO * np.sin(th) * np.sin(PHI)
    h = np.abs(X * np.sin(th) - RHO * np.sin(PHI) * np.cos(th))
    Om = Om(l, h, Rs)
    

    # 被積分関数
    f = 1 / l**2

    # 円柱座標の体積要素
    integrand = f * RHO

    # phi積分
    I_phi = np.trapezoid(integrand, phi, axis=2)

    # rho積分
    I_rho = np.trapezoid(I_phi, rho, axis=1)

    # x積分
    I = np.trapezoid(I_rho, x, axis=0)

    result.append(I)

result = np.array(result)

plt.plot(thr, result)
plt.show()


#nl = L + x[i] * np.cos(thr) + rho[j] * np.sin(thr) * np.sin(phi[k])
#h = np.abs(x[i] * np.sin(thr) - rho[j] * np.sin(phi[k]) * np.cos(phi[k]))
#xhp = x[i] + L * np.cos(thr)
#yhp = rho[j] * np.sin(phi[k]) + L * np.sin(thr)
#zhp = rho[j] * np.cos(phi[k])
#HP = np.sqrt(xhp**2 + yhp**2 + zhp**2)
#yT = yhp * (c + L * np.cos(thr)) / xhp - L * np.sin(thr)
#zT = zhp * (c + L * np.cos(thr)) / xhp