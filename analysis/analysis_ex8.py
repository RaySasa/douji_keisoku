# 同時測定の角度分布(円柱) fitting + cos(theta)(1-alpha sin^2(theta))補正

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ellipk

plt.rcParams["text.usetex"] = True

# 実験データ
xex = np.array([0, 10, 20, 25, 30, 45])
yex = np.array([121.602, 96.827, 33.349, 8.325, 1.804, 1.11])

# 定数
M = 100

L = 5.41
Rs = 1.27
Rb = 2.54
c = 2.7
d = 0.2
r = 0.4


def Omega(l, h, R):
    eta2 = 4 * R * h / ((h + R)**2 + l**2)
    return 2 * np.pi * l / np.sqrt(h**2 + l**2) - 4 * l * ellipk(eta2) / np.sqrt((h + R)**2 + l**2)


def safe_arcsin(u):
    return np.arcsin(np.clip(u, -1, 1))


def Fc(x):
    return x * np.sqrt(np.abs(Rb**2 - x**2)) + Rb**2 * safe_arcsin(x / Rb)


def Fe(x, xT, a, b):
    return b / a * ((x - xT) * np.sqrt(np.abs(a**2 - (x - xT)**2)) + a**2 * safe_arcsin((x - xT) / a))


def calc_distribution(theta_deg):
    theta_deg = np.asarray(theta_deg)
    theta_rad = np.radians(theta_deg)

    x_arr = np.linspace(-d / 2, d / 2, M)
    rho_arr = np.linspace(0, r, M)
    phi_arr = np.linspace(0, 2 * np.pi, M)

    x, rho, phi = np.meshgrid(x_arr, rho_arr, phi_arr, indexing="ij")

    result = []

    for th in theta_rad:
        l = L + x * np.cos(th) + rho * np.sin(th) * np.sin(phi)
        h = np.abs(x * np.sin(th) - rho * np.sin(phi) * np.cos(th))
        Om = Omega(l, h, Rs)

        xhp = x + L * np.cos(th)
        yhp = rho * np.sin(phi) + L * np.sin(th)
        zhp = rho * np.cos(phi)

        yT = yhp * (c + L * np.cos(th)) / xhp - L * np.sin(th)
        zT = zhp * (c + L * np.cos(th)) / xhp
        xT = np.sqrt(yT**2 + zT**2)

        cpsi = xhp / np.sqrt(xhp**2 + yhp**2 + zhp**2)

        b = (c - x) * np.sqrt(xhp**2 + yhp**2 + zhp**2) / xhp * np.sqrt(Om / np.pi)
        a = b / cpsi

        D = (b * xT)**2 - (b**2 - a**2) * (Rb**2 - b**2)
        den = b**2 - a**2

        mask = (D >= 0) & (np.abs(den) > 1e-12)

        x0 = np.full_like(D, np.nan)
        x0[mask] = (b[mask]**2 * xT[mask] - a[mask] * np.sqrt(D[mask])) / den[mask]

        n = np.zeros_like(x)

        mask1 = xT >= Rb + a
        mask2 = xT <= Rb - a
        mask3 = (~mask1) & (~mask2)

        n[mask1] = 0
        n[mask2] = Om[mask2]

        n[mask3] = (Fc(Rb) - Fc(x0[mask3]) + Fe(x0[mask3], xT[mask3], a[mask3], b[mask3]) - Fe(xT[mask3] - a[mask3], xT[mask3], a[mask3], b[mask3])) * cpsi[mask3] * Om[mask3] / (np.pi * b[mask3]**2)

        integrand = n * rho

        I_phi = np.trapezoid(integrand, phi_arr, axis=2)
        I_rho = np.trapezoid(I_phi, rho_arr, axis=1)
        I = np.trapezoid(I_rho, x_arr, axis=0)

        result.append(I)

    return np.array(result)


def correction(theta_deg, alpha):
    theta = np.radians(theta_deg)
    return np.cos(theta) * (1 - alpha * np.sin(theta)**2)


def fit_func(theta_deg, N0, alpha):
    model = calc_distribution(theta_deg)
    return N0 * model / model[0] * correction(theta_deg, alpha)


popt, pcov = curve_fit(fit_func, xex, yex, p0=[yex[0], 0.3], bounds=([0, 0], [np.inf, 1]))

N0_fit, alpha_fit = popt
N0_err, alpha_err = np.sqrt(np.diag(pcov))

print(f"N0 = {N0_fit:.6g} ± {N0_err:.6g}")
print(f"alpha = {alpha_fit:.6g} ± {alpha_err:.6g}")


theta_plot = np.linspace(0, 75, 100)

model_plot = calc_distribution(theta_plot)
y_theory = N0_fit * model_plot / model_plot[0]
y_fit = fit_func(theta_plot, N0_fit, alpha_fit)

plt.figure()
plt.scatter(xex, yex, label="Experiment")
plt.plot(theta_plot, y_theory, linestyle="--", label="Original theory")
plt.plot(theta_plot, y_fit, label=r"$\cos\theta(1-\alpha\sin^2\theta)$ corrected fit")

plt.xlabel(r"$\theta$ [degree]")
plt.ylabel(r"$N(\theta)$")
plt.grid(linestyle="--", linewidth=0.5)
plt.legend()

plt.savefig("tex/analysis_ex10_fit_cos_sin2.pdf", dpi=300, bbox_inches="tight")
plt.show()