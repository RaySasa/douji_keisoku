# 同時測定の角度分布(円柱) fitting

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ellipk

M = 40

L = 5
Rs = 1.27
Rb = 2.54

xex = np.array([0, 10, 20, 25, 30, 35, 45, 60, 75])
yex = np.array([11.281, 10.708, 4.967, 1.521, 1.410, 0.678, 0.408, 0.520, 0.408])


def Omega(l, h, R):
    eta2 = 4 * R * h / ((h + R)**2 + l**2)
    return (
        2*np.pi*l / np.sqrt(h**2 + l**2)
        - 4*l*ellipk(eta2) / np.sqrt((h + R)**2 + l**2)
    )


def safe_arcsin(u):
    return np.arcsin(np.clip(u, -1, 1))


def calc_distribution(theta_deg, c, d, r):
    thr = np.radians(theta_deg)

    x_arr = np.linspace(-d/2, d/2, M)
    rho_arr = np.linspace(0, r, M)
    phi_arr = np.linspace(0, 2*np.pi, M)

    x, rho, phi = np.meshgrid(x_arr, rho_arr, phi_arr, indexing="ij")

    def Fc(u):
        return (
            u*np.sqrt(np.maximum(Rb**2 - u**2, 0))
            + Rb**2 * safe_arcsin(u / Rb)
        )

    def Fe(u, xT, a, b):
        s = u - xT
        return b / a * (
            s*np.sqrt(np.maximum(a**2 - s**2, 0))
            + a**2 * safe_arcsin(s / a)
        )

    result = []

    for th in thr:
        l = L + x*np.cos(th) + rho*np.sin(th)*np.sin(phi)

        h = np.abs(
            x*np.sin(th)
            - rho*np.sin(phi)*np.cos(th)
        )

        Om = Omega(l, h, Rs)

        xT_inside = (
            (
                rho**2
                + L**2*np.sin(th)**2
                + 2*rho*L*np.sin(phi)*np.sin(th)
            )
            * (c + L*np.cos(th))**2
            / (x + L*np.cos(th))**2
            + (L*np.sin(th))**2
            - 2
            * (rho*np.sin(phi) + L*np.sin(th))
            * (c + L*np.cos(th))
            * L*np.sin(th)
            / (x + L*np.cos(th))
        )

        xT = np.sqrt(np.maximum(xT_inside, 0))

        cpsi = (
            (x + L*np.cos(th))
            / np.sqrt(
                x**2 + L**2 + rho**2
                + 2*x*L*np.cos(th)
                + 2*rho*L*np.sin(th)*np.sin(phi)
            )
        )

        b = (
            (c - x)
            * np.sqrt(
                (
                    x**2 + L**2 + rho**2
                    + 2*x*L*np.cos(th)
                    + 2*rho*L*np.sin(th)*np.sin(phi)
                )
                * Om / (2*np.pi)
            )
            / (x + L*np.cos(th))
        )

        a = b / cpsi

        inside_x0 = (
            (b*xT)**2
            - (b**2 - a**2) * (Rb**2 - b**2))

        denom_x0 = b**2 - a**2

        valid_x0 = (
            (inside_x0 >= 0)
            & (np.abs(denom_x0) > 1e-12))

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


def fit_func(theta_deg, A, B, c, d, r):
    theory = calc_distribution(theta_deg, c, d, r)
    return A * theory + B

#def fit_func(theta_deg, A, B, c):
#    theory = calc_distribution(theta_deg, c, d=0.2, r=0.5)
#    return A * theory + B


p0 = [1.0, 0.0, 2.0, 0.2, 0.5]

bounds = (
    [0, -np.inf, 0.0, 0.01, 0.01],
   [np.inf, np.inf, 10.0, 2.0, 2.0]
)

#p0 = [1.0, 0.0, 2.0]

#bounds = (
#    [0, -np.inf, 0.0],
#    [np.inf, np.inf, 20.0]
#)

popt, pcov = curve_fit(
    fit_func,
    xex,
    yex,
    p0=p0,
    bounds=bounds,
    maxfev=20000
)

A_fit, B_fit, c_fit, d_fit, r_fit = popt
#A_fit, B_fit, c_fit = popt
perr = np.sqrt(np.diag(pcov))

print("===== fitting result =====")
print("A =", A_fit, "+/-", perr[0])
print("B =", B_fit, "+/-", perr[1])
print("c =", c_fit, "+/-", perr[2])
print("d =", d_fit, "+/-", perr[3])
print("r =", r_fit, "+/-", perr[4])

theta_plot = np.linspace(0, 75, 300)
fit_y = fit_func(theta_plot, *popt)

plt.scatter(xex, yex, label="experiment")
plt.plot(theta_plot, fit_y, label="fit")
plt.xlabel("theta [deg]")
plt.ylabel("counts")
plt.legend()
plt.savefig("tex/analysis_ex2_fit.pdf", dpi=300, bbox_inches="tight")
plt.show()