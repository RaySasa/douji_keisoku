#課題4 FWHMとエネルギー分解能

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#texフォントを使用
plt.rcParams["text.usetex"] = True

# 線源名
labels = [
    r"$^{22}\mathrm{Na}$",
    r"$^{22}\mathrm{Na}$",
    r"$^{57}\mathrm{Co}$",
    r"$^{60}\mathrm{Co}$",
    r"$^{60}\mathrm{Co}$",
    r"$^{137}\mathrm{Cs}$"
]

#値の引き継ぎ
x_p = np.array([336, 138, 35, 308, 350])
y_p = np.array([1275, 511, 122, 1173, 1333])
(p_p, cov_p) = np.polyfit(x_p, y_p, 1, cov=True)
a_p = p_p[0]
b_p = p_p[1]
da_p = np.sqrt(cov_p[0, 0])
db_p = np.sqrt(cov_p[1, 1])

#data
x_ch = np.array([13.5, 21.5, 4.5, 16.5, 20.0, 16.0])#FWHM(ch)
y_ch = np.array([138, 336, 35, 308, 350, 178])#peak(ch)

#keVに変換
x_keV = x_ch * a_p
y_keV = y_ch * a_p + b_p

#分解能Rを計算
R = x_keV / y_keV

#\frac{1}{\sqrt{E}}
Y = 1 / np.sqrt(y_keV)

#最小2乗法(直線)
(p_l, cov_l) = np.polyfit(Y, R, 1, cov=True)
a_l = p_l[0]
b_l = p_l[1]
da_l = np.sqrt(cov_l[0, 0])
db_l = np.sqrt(cov_l[1, 1])
R_lfit = a_l * Y + b_l

#最小2乗法(曲線)
#フィット関数
def f(x, a, b):
    return a /(np.sqrt(x)) + b
#フィット
popt, pcov = curve_fit(f, y_keV, R)
#係数
a_c = popt[0]
b_c = popt[1]
#誤差
da_c = np.sqrt(pcov[0,0])
db_c = np.sqrt(pcov[1,1])
#曲線
y_cfit = np.linspace(min(y_keV), max(y_keV), 300)
R_cfit = f(y_cfit, a_c, b_c)


#描画
fig, ax = plt.subplots(1, 2, figsize = (12,5))

ax[0].scatter(y_keV, R, label="data", marker='o', s=20)
ax[0].plot(y_cfit, R_cfit, label="fit")

ax[1].scatter(Y, R, label="data", marker='o', s=20)
ax[1].plot(Y, R_lfit)

#ラベル
for i in range(0,2):
    ax[i].legend()
    ax[i].grid(linestyle = "--", linewidth = 0.5)

#画像を保存
plt.savefig("tex/analysis_1.pdf", dpi=300, bbox_inches="tight")
plt.show()