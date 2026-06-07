#TACの時間較正 11)

import matplotlib.pyplot as plt
import numpy as np

#texフォントを使用
plt.rcParams["text.usetex"] = False

x = np.array([25, 50, 76, 102, 128, 154, 180, 206, 232, 259])
y = np.arange(56, 507, 50)

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
    rf"$y = ({a:.3f} \pm {da:.3f})x + ({b:.1f} \pm {db:.1f})$"
)

#点をプロット
plt.scatter(x, y, label="data", marker='o', s=20)

#直線
plt.plot(x, y_fit,label=fit_label, linewidth = 1.0)

#ラベル
plt.xlabel("delay time (ch)")
plt.ylabel("delay time (ns)")
plt.legend()
plt.grid(linestyle = "--", linewidth = 0.5)

#画像を保存
plt.savefig("tex/analysis_3.pdf", dpi=300, bbox_inches="tight")
plt.show()