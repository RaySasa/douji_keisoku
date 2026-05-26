#コンプトン散乱の散乱断面積の分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#データ
xex = np.array([30, 45, 60, 75, 90])#角度
yex = np.array([3.01, 3.208, 2.868, 2.431, 2.072])#cps

r = np.array([0.015, 0.02, 0.03, 0.04, 0.07])
x = 2.54
rho = 3.67
mu = r * 3.67

ycor = yex / (1 - np.exp(- mu * x))

# constants [keV]
E_gamma = 661.7
mc2 = 511.0

alpha = E_gamma / mc2

def klein_nishina(theta, a):
    """
    theta: scattering angle [rad]
    returns y = dσ/dΩ up to factor r_e^2
    """
    return a * (
        (1 / (1 + alpha * (1 - np.cos(theta))))**2
        * ((1 + np.cos(theta)**2) / 2)
        * (1 + (alpha**2 * (1 - np.cos(theta))**2)
           / ((1 + np.cos(theta)**2) * (1 + alpha * (1 - np.cos(theta)))))
    ) 

# theta [degree]
theta_deg = np.linspace(0, 90, 1000)
theta_rad = np.deg2rad(theta_deg)

y = klein_nishina(theta_rad, 1)

#フィット
popt, pcov = curve_fit(klein_nishina, xex, ycor,  p0=[40],bounds=(40, np.inf))
#係数
a = popt[0]
#b = popt[1]
#誤差
da = np.sqrt(pcov[0,0])
#db = np.sqrt(pcov[1,1])
#曲線
y_fit = klein_nishina(theta_rad, a)

print(ycor)

plt.plot(theta_deg, y_fit)
plt.xlabel("theta [degree]")#
plt.ylabel("y")
plt.grid()

#plt.scatter(xex, yex, label="data", marker='o', s=20)
plt.scatter(xex, ycor, label="data", marker='o', s=20)
plt.savefig("tex/analysis_6.pdf", dpi=300, bbox_inches="tight")
plt.show()