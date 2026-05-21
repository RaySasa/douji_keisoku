#同時測定の角度分布(点近似)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def Om(l,Rs):
    return 2 * np.pi * (1 - l / np.sqrt(Rs**2 + l**2))

th = np.linspace(0,  80, 100000)
thr = np.radians(th)

l = 5
Rs = 1.27
Rb = 2.54
c = 3
Om = Om(l, Rs)
b = c * np.sqrt(Om / (2 * np.pi)) / np.cos(thr)
a = b / np.cos(thr)
xT = c * np.tan(thr)
x0 = np.sqrt(np.abs(c**2 + Rb**2 - b**2)) / np.sin(thr) - c / np.tan(thr)

def Fc(x):
    return x*np.sqrt(np.abs(Rb**2 - x**2)) + Rb**2 * np.arcsin(x / Rb)

def Fe(i, x):
    return b[i] / a[i] * ( (x - xT[i]) * np.sqrt(np.abs(a[i]**2 - (x - xT[i])**2)) + a[i]**2 * np.arcsin((x - xT[i]) / a[i]) )

S = np.zeros(len(thr))

for i in range(len(thr)):
    if xT[i] >= Rb +a[i]:
        S[i] = 0
    elif xT[i] <= Rb - a[i]:
        S[i] = np.pi * a[i] * b[i]
    elif xT[i] <= Rb +a[i] and xT[i] >= Rb - a[i]:
        S[i] = Fc(Rb) - Fc(x0[i]) + Fe(i, x0[i]) - Fe(i, xT[i] - a[i])
    #else:
        #S[i] = np.pi * b[i]**2 / np.cos(thr[i])
F = S * np.cos(thr) / (np.pi * b**2)
M = np.zeros(len(thr))
for i in range(len(thr)):
    M = Fc(Rb) - Fc(x0[i]) + Fe(i, x0[i]) - Fe(i, xT[i] - a[i])

#plt.plot(th, F)
#plt.plot([0,80],[0,0])    
plt.plot(th, np.full([len(th), 1], Rb) , label= "R")
plt.plot(th, xT + a, label= "x+a")
plt.plot(th, xT - a, label= "x-a")
plt.legend()
plt.savefig("tex/analysis_ex.pdf", dpi=300, bbox_inches="tight")
plt.show()