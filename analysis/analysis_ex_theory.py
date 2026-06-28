#同時測定の角度分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import *

#texフォントを使用
plt.rcParams["text.usetex"] = False

#変数の設定
M = 100
the = np.linspace(0, 80, M)
thr = np.radians(the)

L = 10000
c = 10000

#固定パラメータ
Rs = 1.27
Rb = 2.54
d = 0.2
r = 0.4

#積分変数の設定
x_arr = np.linspace(-d/2, d/2, M)
rho_arr = np.linspace(0, r, M)
phi_arr = np.linspace(0, 2*np.pi, M)
x, rho, phi = np.meshgrid(x_arr, rho_arr, phi_arr, indexing="ij")

#立体角の定義
def Omega(l, h, R):
    eta2 = 4 * R * h/((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt( (h + R)**2 + l**2 )

#arcsinの計算の安定化
def safe_arcsin(u):
    return np.arcsin(np.clip(u, -1, 1))

#原始関数
def Fc(x):
    return x*np.sqrt(np.abs(Rb**2 - x**2)) + Rb**2 * safe_arcsin(x / Rb)

def Fe(x, xT, a, b):
    return b / a * ( (x - xT) * np.sqrt(np.abs(a**2 - (x - xT)**2)) + a**2 * safe_arcsin((x - xT) / a) )

# 結果保存
result = []

#それぞれのthetaに対して計算
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

    cpsi = xhp / np.sqrt(xhp**2 + yhp**2 + zhp**2)
    b = (c - x) * np.sqrt(xhp**2 + yhp**2 + zhp**2) /xhp * np.sqrt(Om / np.pi)
    a = b / cpsi 

    D = (b * xT)**2 - (b**2 - a**2) * (Rb**2 - b**2)
    den = b**2 - a**2

    mask = (D >= 0) & (np.abs(den) > 1e-12)
    

    x0 = np.full_like(D, np.nan)

    x0[mask] = (b[mask]**2 * xT[mask] - a[mask] * np.sqrt(D[mask])) / den[mask]


    #被積分関数
    n = np.zeros_like(x)
    
    #場合わけ
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

#検出数N
N = np.array(result)


#グラフを表示
plt.xlabel(rf"$\theta$ [degree]")
plt.ylabel(rf"$N(\theta)/N_0$")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.plot(the, N/ N[0])
plt.savefig("tex/analysis_ex_L=c.pdf", dpi=300, bbox_inches="tight")
plt.show()