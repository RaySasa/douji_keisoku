#コンプトン散乱のエネルギー分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#エネルギー較正
x0 = np.array([60, 236])
y0 = np.array([122, 511])

a = (y0[1]-y0[0])/(x0[1]-x0[0])
b = y0[0] - a * x0[0]

#chのデータ
x = np.array([30, 45, 60, 75, 90])
ych = np.array([256, 219, 187, 158, 136])

ykeV = ych * a + b

#エラーバー
count = np.array([287.665, 962.366, 860.544, 729.383, 621.670])
sigma = np.array([31.957, 62.833, 63.785, 44.105, 46.461]) / 2 * a
y_err = sigma / np.sqrt(count)

#エネルギー分布の理論値
e = 661.7 #入射γ線のエネルギー
mc2 = 511 #定数

def E(t):
    return e / ( 1 + e * ( 1 - np.cos(t * 2 * np.pi / 360) )/ mc2)

x_fit = np.linspace(0, 90, 300)
Ec = E(x_fit)

print(ych,ykeV,y_err)

#plt.scatter(x, ykeV, label="data", marker='o', s=20)
plt.errorbar(x, ykeV, yerr = y_err, capsize=5, fmt='.', markersize=0, label= 'data')
plt.plot(x_fit, Ec, linewidth = 1.0, label= 'theory')


plt.legend()
plt.xlabel(r"$\theta\,[degree]$")
plt.ylabel("Gamma-ray Energy")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.savefig("tex/analysis_5.pdf", dpi=300, bbox_inches="tight")
plt.show()