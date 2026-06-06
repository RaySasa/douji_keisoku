#課題4 FWHMとエネルギー分解能

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#texフォントを使用
plt.rcParams["text.usetex"] = False

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

#誤差の計算
dx = x_ch * da_p
dy = np.sqrt(y_ch**2 * da_p**2 + db_p**2)

#分解能Rを計算
R = x_keV / y_keV
#誤差
dR = np.sqrt( (dx/y_keV)**2 + (dy*x_keV/y_keV**2)**2) 
print(x_ch)
print(R)
print(dR)

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
#popt, pcov = curve_fit(f, y_keV, R)
#係数
#a_c = popt[0]
#b_c = popt[1]
#誤差
#da_c = np.sqrt(pcov[0,0])
#db_c = np.sqrt(pcov[1,1])
#曲線
y_cfit = np.linspace(min(y_keV), max(y_keV), 300)
R_cfit = f(y_cfit, a_l, b_l)

#描画
fig, ax = plt.subplots(1, 2, figsize = (12,5))

#曲線のラベル
fit_label_0 = (
    rf"$y = ({a_l:.3f} \pm {da_l:.3f})/\sqrt{{x}} + ({b_l:.3f} \pm {db_l:.3f})$"
)
fit_label_1 = (
    rf"$y = ({a_l:.3f} \pm {da_l:.3f})x + ({b_l:.3f} \pm {db_l:.3f})$"
)

# 各点にラベル
for i in range(len(R)):
    if i == 2:
        ax[0].annotate(
            labels[i],
            (y_keV[i], R[i]),
            textcoords="offset points",
         xytext=(5,0),
            fontsize=12
        )
    else:    
        ax[0].annotate(
            labels[i],
            (y_keV[i], R[i]),
            textcoords="offset points",
         xytext=(-30,-10),
            fontsize=12
        )

for i in range(len(R)):
    if i == 3:
        ax[1].annotate(
            labels[i],
            (Y[i], R[i]),
            textcoords="offset points",
        xytext=(5,-8),
            fontsize=12
        )
    elif i == 2:
        ax[1].annotate(
            labels[i],
            (Y[i], R[i]),
            textcoords="offset points",
         xytext=(-10,-20),
            fontsize=12
        )
    else:
        ax[1].annotate(
            labels[i],
            (Y[i], R[i]),
            textcoords="offset points",
         xytext=(5,-5),
            fontsize=12
        )
               
#点と曲線
ax[0].scatter(y_keV, R, label="data", marker='o', s=20)
ax[0].plot(y_cfit, R_cfit, label=fit_label_0, linewidth = 1.0)

ax[1].scatter(Y, R, label="data", marker='o', s=20)
ax[1].plot(Y, R_lfit,label=fit_label_1, linewidth = 1.0)

#ラベル
ax[0].set_xlabel("E [keV]")
ax[1].set_xlabel("E^(-1/2) [keV^(-1/2)]")

for i in range(0,2):
    ax[i].legend()
    ax[i].grid(linestyle = "--", linewidth = 0.5)
    ax[i].set_ylabel("Resolution")

#画像を保存
plt.savefig("tex/analysis_2.pdf", dpi=300, bbox_inches="tight")
plt.show()
