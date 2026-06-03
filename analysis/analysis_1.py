#課題1 チャンネルとγ線のエネルギーの関係

import matplotlib.pyplot as plt
import numpy as np

#texフォントを使用
plt.rcParams["text.usetex"] = True


#データ
x = np.array([336, 138, 35, 308, 350])#チャンネル
y = np.array([1275, 511, 122, 1173, 1333])#エネルギー

# 線源名
labels = [
    r"$^{22}\mathrm{Na}$",
    r"$^{22}\mathrm{Na}$",
    r"$^{57}\mathrm{Co}$",
    r"$^{60}\mathrm{Co}$",
    r"$^{60}\mathrm{Co}$"
]


#最小2乗法
(p, cov) = np.polyfit(x, y, 1, cov=True)

#係数
a = p[0]
b = p[1]

#誤差
da = np.sqrt(cov[0, 0])
db = np.sqrt(cov[1, 1])

#フィッティング
y_fit = a * x + b

#描画
fit_label = (
    rf"$y = ({a:.2f} \pm {da:.2f})x + ({b:.0f} \pm {db:.0f})$"
)

#Csの点
x1 = 178
y1 = 661.7

#点をプロット
plt.scatter(x, y, label="data", marker='o', s=20)
#plt.scatter(178, 661.7, marker='o', s=20, color='tab:blue')

# 各点にラベル
for i in range(len(x)):
    if i == 2:
        plt.annotate(
            labels[i],
            (x[i], y[i]),
            textcoords="offset points",
         xytext=(-15,10),
            fontsize=12
        )
    else:    
        plt.annotate(
            labels[i],
            (x[i], y[i]),
            textcoords="offset points",
         xytext=(-30,0),
            fontsize=12
        )

#Csのラベル
#plt.annotate(r"$^{137}\mathrm{Cs}$", (178, 661.7), textcoords="offset points", xytext=(-30,0), fontsize=12)

#直線
plt.plot(x, y_fit, label=fit_label, linewidth = 1.0)

#ラベル
plt.xlabel("Gamma-ray Energy [ch]}")
plt.ylabel("Gamma-ray Energy [keV]")
plt.legend()
plt.grid(linestyle = "--", linewidth = 0.5)

#画像を保存
plt.savefig("tex/analysis_1.pdf", dpi=300, bbox_inches="tight")
plt.show()