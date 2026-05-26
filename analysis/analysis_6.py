#コンプトン散乱の散乱断面積の分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#データ
xex = np.array([45, 60, 75, 90])#角度
yex = np.array([3.208, 2.868, 2.431, 2.072])#cps

r = np.array([0.085, 0.1, 0.105, 0.16])
x = 2.54
rho = 3.67
mu = r * 3.67

ycor = yex * (1 - np.exp(- mu * x))

# constants [keV]
E_gamma = 661.7
mc2 = 511.0

alpha = E_gamma / mc2

def klein_nishina(theta):
    """
    theta: scattering angle [rad]
    returns y = dσ/dΩ up to factor r_e^2
    """
    return (
        (1 / (1 + alpha * (1 - np.cos(theta))))**2
        * ((1 + np.cos(theta)**2) / 2)
        * (1 + (alpha**2 * (1 - np.cos(theta))**2)
           / ((1 + np.cos(theta)**2) * (1 + alpha * (1 - np.cos(theta)))))
    )

# theta [degree]
theta_deg = np.linspace(0, 90, 1000)
theta_rad = np.deg2rad(theta_deg)

y = klein_nishina(theta_rad)

def fit_func(theta_rad, a, b):
    return a * klein_nishina(theta_rad) + b

#フィット
popt, pcov = curve_fit(fit_func, a, b)
#係数
#a_c = popt[0]
#b_c = popt[1]
#誤差
#da_c = np.sqrt(pcov[0,0])
#db_c = np.sqrt(pcov[1,1])
#曲線
y_cfit = np.linspace(min(y_keV), max(y_keV), 300)
R_cfit = f(y_cfit, a_l, b_l)

#plt.plot(theta_deg, y)
plt.xlabel("theta [degree]")
plt.ylabel("y")
plt.grid()

plt.scatter(xex, yex, label="data", marker='o', s=20)
#plt.scatter(xex, ycor, label="data", marker='o', s=20)

plt.show()