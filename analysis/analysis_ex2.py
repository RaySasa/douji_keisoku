#同時測定の角度分布(円柱)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import *

#texフォントを使用
plt.rcParams["text.usetex"] = True

M = 100
the = np.linspace(0, 90, 100)
thr = np.radians(the)

L = 5.5
Rs = 1.27
Rb = 2.54
c = 5.5
d = 0.5
r = 0.5
lam = 0.7 #減衰率

x_arr = np.linspace(-d/2, d/2, M)
rho_arr = np.linspace(0, r, M)
phi_arr = np.linspace(0, 2*np.pi, M)
x, rho, phi = np.meshgrid(x_arr, rho_arr, phi_arr, indexing="ij")

def Omega(l, h, R):
    eta2 = 4 * R * h/((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt( (h + R)**2 + l**2 )

def safe_arcsin(u):
    return np.arcsin(np.clip(u, -1, 1))

def Fc(x):
    return x*np.sqrt(np.abs(Rb**2 - x**2)) + Rb**2 * safe_arcsin(x / Rb)

def Fe(x, xT, a, b):
    return b / a * ( (x - xT) * np.sqrt(np.abs(a**2 - (x - xT)**2)) + a**2 * safe_arcsin((x - xT) / a) )

# 結果保存
result = []

for th in thr:

    # 色々定義
    l = L + x * np.cos(th) + rho * np.sin(th) * np.sin(phi)
    h = np.abs(x * np.sin(th) - rho * np.sin(phi) * np.cos(th))
    Om = Omega(l, h, Rs)

    xhp = x + L*np.cos(th)
    yhp = rho * np.sin(phi) + L * np.sin(th)
    zhp = rho * np.cos(phi)
    yT = yhp * (c + L * np.cos(th))/xhp - L * np.sin(th)
    zT = zhp * (c + L * np.cos(th))/xhp
    xT = np.sqrt(yT**2 + zT**2)
    #xT = np.sqrt( (rho**2 + L**2 * np.sin(th)**2 + 2 * rho * L * np.sin(phi) * np.sin(th))*(c + L * np.cos(th))**2 / (x + L * np.cos(th))**2 + (L * np.sin(th))**2 - 2 * (rho * np.sin(phi) + L * np.sin(th)) * (c + L * np.cos(th)) * L * np.sin(th) / (x + L * np.cos(th))) 
    cpsi = xhp / np.sqrt(xhp**2 + yhp**2 + zhp**2)
    #cpsi = (x + L * np.cos(th)) / (np.sqrt( (x**2 + L**2 +rho**2 + 2 * x * L * np.cos(th) + 2 * rho * L * np.sin(th) * np.sin(phi))))
    b = (c - x) * np.sqrt(xhp**2 + yhp**2 + zhp**2) /xhp * np.sqrt(Om / np.pi)
    #b = (c - x) * np.sqrt( (x**2 + L**2 +rho**2 +2 * x * L * np.cos(th) + 2 * rho * L * np.sin(th) * np.sin(phi)) * Om / (np.pi)) / (x + L * np.cos(th))
    a = b / cpsi 
    #a = (x**2 + L**2 +rho**2 +2 * x * L * np.cos(th) + 2 * rho * L * np.sin(th) * np.sin(phi)) * (c - x) * Om / (2 * np.pi * (x + L * np.cos(th))**2)
    x0 = (b**2 * xT - a * np.sqrt((b * xT)**2 - (b**2 - a**2) * (Rb**2 - b**2))) / (b**2 - a**2)

    #被積分関数
    n = np.zeros_like(x)
    
    mask1 = xT >= Rb + a
    mask2 = xT <= Rb - a
    mask3 = (~mask1) & (~mask2)

    n[mask1] = 0
    n[mask2] = Om[mask2]
    n[mask3] = (Fc(Rb) - Fc(x0[mask3]) + Fe(x0[mask3], xT[mask3], a[mask3], b[mask3]) - Fe(xT[mask3] - a[mask3], xT[mask3], a[mask3], b[mask3])) * cpsi[mask3] * Om[mask3] / (np.pi * b[mask3]**2) 


    # 円柱座標の体積要素
    integrand = n * rho

    # phi積分
    I_phi = np.trapezoid(integrand, phi_arr, axis=2)

    # rho積分
    I_rho = np.trapezoid(I_phi, rho_arr, axis=1)

    # x積分
    I = np.trapezoid(I_rho, x_arr, axis=0)

    result.append(I)

result = np.array(result)


plt.xlabel(rf"$\theta$ [degree]")
plt.ylabel(rf"$N(\theta)/N_0$")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.plot(the, result/ result[0])
# 枠線を表示
#ax = plt.gca()
#ax.spines['top'].set_visible(True)
#ax.spines['right'].set_visible(True)
#ax.spines['bottom'].set_visible(True)
#ax.spines['left'].set_visible(True)
plt.savefig("tex/analysis_ex2_final.pdf", dpi=300, bbox_inches="tight")
plt.show()