#問2 角度分布

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

T = 300
x = np.array([0, 10, 20, 25, 30, 35, 45, 60, 75])
y = np.array([11.281, 10.708, 4.967, 1.521, 1.410, 0.678, 0.408, 0.520, 0.408])

x1 = np.array([0, 5, 10, 15, 20, 25, 30, 35, 45, 60])
y1 = np.array([145.691, 128.762, 94.860, 78.117, 41.344, 15.050, 2.908, 1.127, 0.762, 0.742])

yer = np.array([])
for i in range(len(y)):
    if i == 2:
        yer = np.append(yer, np.sqrt(y[i] * 280) / 280)
    elif i == 3:
        yer = np.append(yer, np.sqrt(y[i] *243) / 243)
    else:
        yer = np.append(yer, np.sqrt(y[i] * T) / T)

plt.scatter(x1, y1, label="data", marker='o', s=20)
plt.savefig("tex/analysis_4.pdf", dpi=300, bbox_inches="tight")
plt.show()