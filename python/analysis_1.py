#課題1 チャンネルとγ線のエネルギーの関係

import matplotlib.pyplot as plt
import numpy as np


#データ
x = np.array([336, 138, 35, 308, 350])#チャンネル
y = np.array([1275, 511, 122, 1173, 1333])#エネルギー

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
plt.scatter(x, y, label="data")
plt.plot(x, y_fit, label="fit")

plt.xlabel("ch")
plt.ylabel("keV")
plt.legend()
plt.grid(linestyle = "--", linewidth =0.5)

plt.show()