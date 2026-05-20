#同時測定の角度分布(点近似)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def Om(l,Rs):
    return 2 * np.pi * (1 - l / np.sqrt(Rs**2 + l**2))





print(Om(5,2.27))


