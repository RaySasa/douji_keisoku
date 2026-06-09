#実験3のエネルギー較正

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# 線源名
labels = [
    r"$^{22}\mathrm{Na}$",
    r"$^{57}\mathrm{Co}$",
    r"$^{137}\mathrm{Cs}$"
]


#エネルギー較正
x0 = np.array([60, 236])
y0 = np.array([122, 511])

x1 = 300
y1 = 661.7

a = (y0[1]-y0[0])/(x0[1]-x0[0])
b = y0[0] - a * x0[0]

x = np.linspace(0, 350)
y = a * x + b

#描画
fit_label = (
    rf"$y = {a:.2f}x {b:.1f}$"
)

x3 = np.append(x0, x1)
y3 = np.append(y0, y1)

# 各点にラベル
for i in range(len(x3)):
    plt.annotate(
        labels[i],
        (x3[i], y3[i]),
        textcoords="offset points",
        xytext=(-30,5),
        fontsize=12
        )

plt.scatter(x0, y0, label= 'data', s= 20)
plt.scatter(x1, y1, s= 20)
plt.plot(x, y, label= fit_label)

plt.xlabel("Gamma-ray Energy (ch)")
plt.ylabel("Gamma-ray Energy (keV)")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.legend()
plt.savefig("tex/analysis_9.pdf", dpi=300, bbox_inches="tight")
plt.show()