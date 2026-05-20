#コンプトン散乱の散乱断面積の分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#データ
x = np.array([30, 45, 60, 75, 90])#角度
y = np.array([0.959, 3.208, 2.868, 2.431, 2.072])#cps

plt.scatter(x, y, label="data", marker='o', s=20)
plt.show()