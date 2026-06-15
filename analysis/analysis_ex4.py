# 同時測定の角度分布(円柱) AとBだけfitting

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ellipk

#texフォントを使用
plt.rcParams["text.usetex"] = False

M = 50

L = 5.41
Rs = 1.27
Rb = 2.54

# 固定パラメータ
c_fixed = 5.93
d_fixed = 0.2
r_fixed = 0.4

#新しいデータ L = 5.41
xex = np.array([0, 5, 10, 15, 20, 25, 30, 35, 45, 60])
yex = np.array([145.691, 128.762, 94.860, 78.117, 41.344, 15.050, 2.908, 1.127, 0.762, 0.742])
count = np.array([13.112e3, 10.687e3, 5.692e3, 5.468e3, 5.044e3, 3.281e3, 424.535, 203.925, 189.640, 222.549])
dx = np.full(len(xex), 0.5)
#dx = np.full(len(xex), np.degrees(np.arctan(Rs / L)))
dy = np.sqrt(count) * yex /count

# L = 10.41
#xex = np.array([0, 5, 10, 15, 20, 25, 30, 35, 45])
#yex = np.array([43.539, 43.004, 41.696, 37.384, 24.433, 3.345, 0.519, 0.275, 0.214])
#count = np.array([3.744e3, 6.365e3, 8.172e3, 9.122e3, 2.565e3, 428.115, 110.028, 82.495, 46.628])
#dx = np.full(len(xex), 0.5)
#dy = np.sqrt(count) * yex / count


# 積分用メッシュ
x_arr = np.linspace(-d_fixed / 2, d_fixed / 2, M)
rho_arr = np.linspace(0, r_fixed, M)
phi_arr = np.linspace(0, 2 * np.pi, M)

x, rho, phi = np.meshgrid(
    x_arr,
    rho_arr,
    phi_arr,
    indexing="ij"
)


def safe_arcsin(u):
    return np.arcsin(np.clip(u, -1, 1))


def Omega(l, h, R):
    h = np.maximum(h, 1e-12)

    eta2 = 4 * R * h / ((h + R)**2 + l**2)
    eta2 = np.clip(eta2, 0, 1 - 1e-12)

    return (
        2 * np.pi * l / np.sqrt(h**2 + l**2)
        - 4 * l * ellipk(eta2) / np.sqrt((h + R)**2 + l**2)
    )


def calc_distribution(theta_deg):
    thr = np.radians(theta_deg)
    result = []

    def Fc(u):
        return (
            u * np.sqrt(np.maximum(Rb**2 - u**2, 0))
            + Rb**2 * safe_arcsin(u / Rb)
        )

    def Fe(u, xT, a, b):
        s = u - xT
        return b / a * (
            s * np.sqrt(np.maximum(a**2 - s**2, 0))
            + a**2 * safe_arcsin(s / a)
        )

    for th in thr:
        sin_th = np.sin(th)
        cos_th = np.cos(th)

        denominator = x + L * cos_th
        denominator = np.where(
            np.abs(denominator) < 1e-12,
            1e-12,
            denominator
        )

        distance2 = (
            x**2
            + L**2
            + rho**2
            + 2 * x * L * cos_th
            + 2 * rho * L * sin_th * np.sin(phi)
        )

        distance2 = np.maximum(distance2, 1e-12)

        l = L + x * cos_th + rho * sin_th * np.sin(phi)

        h = np.abs(
            x * sin_th
            - rho * np.sin(phi) * cos_th
        )

        Om = Omega(l, h, Rs)

        xT_inside = (
            (
                rho**2
                + L**2 * sin_th**2
                + 2 * rho * L * np.sin(phi) * sin_th
            )
            * (c_fixed + L * cos_th)**2
            / denominator**2
            + (L * sin_th)**2
            - 2
            * (rho * np.sin(phi) + L * sin_th)
            * (c_fixed + L * cos_th)
            * L * sin_th
            / denominator
        )

        xT = np.sqrt(np.maximum(xT_inside, 0))

        cpsi = denominator / np.sqrt(distance2)

        b = (
            (c_fixed - x)
            * np.sqrt(distance2 * Om / (np.pi))
            / denominator
        )

        b = np.maximum(b, 1e-12)

        a = b / cpsi
        a = np.maximum(a, 1e-12)

        inside_x0 = (
            (b * xT)**2
            - (b**2 - a**2) * (Rb**2 - b**2)
        )

        denom_x0 = b**2 - a**2

        valid_x0 = (
            (inside_x0 >= 0)
            & (np.abs(denom_x0) > 1e-12)
        )

        x0 = np.full_like(x, np.nan)

        x0[valid_x0] = (
            b[valid_x0]**2 * xT[valid_x0]
            - a[valid_x0] * np.sqrt(inside_x0[valid_x0])
        ) / denom_x0[valid_x0]

        n = np.zeros_like(x)

        mask1 = xT >= Rb + a
        mask2 = xT <= Rb - a
        mask3 = (~mask1) & (~mask2) & valid_x0

        n[mask1] = 0
        n[mask2] = Om[mask2]

        n[mask3] = (
            Fc(Rb)
            - Fc(x0[mask3])
            + Fe(x0[mask3], xT[mask3], a[mask3], b[mask3])
            - Fe(
                xT[mask3] - a[mask3],
                xT[mask3],
                a[mask3],
                b[mask3]
            )
        ) * cpsi[mask3] * Om[mask3] / (np.pi * b[mask3]**2)

        integrand = n * rho

        I_phi = np.trapezoid(integrand, phi_arr, axis=2)
        I_rho = np.trapezoid(I_phi, rho_arr, axis=1)
        I = np.trapezoid(I_rho, x_arr, axis=0)

        result.append(I)

    return np.array(result)


def fit_func(theta_deg, A, B):
    theory = calc_distribution(theta_deg)
    return A * theory + B


# 初期値
p0 = [1.0, 0.0]

bounds = (
    [0, 0],
    [np.inf, np.inf]
)

popt, pcov = curve_fit(
    fit_func,
    xex,
    yex,
    p0=p0,
    bounds=bounds,
    maxfev=20000
)

A_fit, B_fit = popt
perr = np.sqrt(np.diag(pcov))
dA, dB = perr

print("===== fitting result =====")
print("A =", A_fit, "+/-", perr[0])
print("B =", B_fit, "+/-", perr[1])
print("c =", c_fixed, "fixed")
print("d =", d_fixed, "fixed")
print("r =", r_fixed, "fixed")

theta_plot = np.linspace(0, 75, 100)
fit_y = fit_func(theta_plot, *popt)

#曲線のラベル
fit_label = (
    rf"cps$ = ({round(A_fit, -2):.0f} \pm {round(dA, -1):.0f})N(\theta)/n"
    rf" + ({round(B_fit, 0):.0f} \pm {round(dB, 0):.0f})$"
)

plt.errorbar(xex, yex, xerr = dx, yerr = dy, capsize=5, fmt='o', label= "data", markersize= 3, color= "crimson", ecolor = "tab:blue")
plt.plot(theta_plot, fit_y, label=fit_label, color= "darkorange")

plt.xlabel(rf"$\theta\,[deg]$")
plt.ylabel("cps [1/s]")
plt.legend()
plt.grid(linestyle = "--", linewidth = 0.5)

plt.savefig("tex/analysis_ex_fit_L5.pdf", dpi=300, bbox_inches="tight")
plt.show()