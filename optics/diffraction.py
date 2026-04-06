"""
Diffraction analysis.
TP Optics L3, Sorbonne.
"""

import numpy as np
from scipy.optimize import curve_fit


def linear_func(x, a, b):
    return a * x + b

def proportional_func(x, a):
    return a * x


def analyze_single_slit(x_pos, orders, D, lam=578e-9):
    """
    Calculate slit width 'a' from minima positions.
    Formula: x = (lam * D / a) * p
    """
    popt, pcov = curve_fit(linear_func, orders, x_pos)

    slope = popt[0]        # = lam * D / a
    intercept = popt[1]
    u_slope = np.sqrt(pcov[0][0])

    a = (lam * D) / slope
    u_a = a * (u_slope / slope)

    residuals = x_pos - linear_func(orders, slope, intercept)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((x_pos - np.mean(x_pos))**2)
    r2 = 1 - ss_res / ss_tot

    return {"a": a, "u_a": u_a, "R2": r2, "slope": slope, "intercept": intercept}


def analyze_double_slit(x_pos, orders, D, lam=578e-9):
    """
    Calculate slit separation 'b' from fringe positions.
    x = (lam * D / b) * n
    """
    popt, pcov = curve_fit(linear_func, orders, x_pos)

    slope = popt[0]
    u_slope = np.sqrt(pcov[0][0])

    b = (lam * D) / slope
    u_b = b * (u_slope / slope)

    return b, u_b


def get_wavelength_grating(angles_deg, orders, d):
    """
    Find wavelength from grating equation: d * sin(theta) = n * lam
    => sin(theta) = (lam / d) * n
    """
    sin_theta = np.sin(np.deg2rad(angles_deg))

    popt, pcov = curve_fit(proportional_func, orders, sin_theta)

    slope = popt[0]        # = lam / d
    u_slope = np.sqrt(pcov[0][0])

    lam = slope * d
    u_lam = u_slope * d

    return lam, u_lam


def get_grating_constant(sin_theta, wavelengths):
    """
    Find grating constant d from known spectral lines (1st order).
    sin(theta) = (1/d) * lam
    Used for mercury lamp in TP5.
    """
    popt, pcov = curve_fit(proportional_func, wavelengths, sin_theta)

    A = popt[0]            # = 1/d
    u_A = np.sqrt(pcov[0][0])

    d = 1 / A
    u_d = u_A / A**2

    return d, u_d


def sinc_intensity_profile(x, a, lam, f, I0):
    """
    Single slit intensity: I = I0 * sinc^2(a*x / (lam*f))
    np.sinc(X) = sin(pi*X) / (pi*X), so X = a*x / (lam*f).
    """
    X = (a * x) / (lam * f)
    return I0 * np.sinc(X)**2


def fringe_visibility(I_max, I_min):
    """V = (Imax - Imin) / (Imax + Imin)."""
    if I_max + I_min == 0:
        return 0.0
    return (I_max - I_min) / (I_max + I_min)


if __name__ == "__main__":
    print("testing diffraction")

    orders = np.array([-3, -2, -1, 1, 2, 3])
    D = 0.3  # 30 cm
    x_pos = np.array([-0.0013, -0.00085, -0.00042, 0.00044, 0.00088, 0.00132])

    res = analyze_single_slit(x_pos, orders, D)
    print(f"a = {res['a']*1e6:.1f} um")
    print(f"u_a = {res['u_a']*1e6:.1f} um")
    print(f"R2 = {res['R2']:.6f}")

    import matplotlib.pyplot as plt
    x = np.linspace(-0.005, 0.005, 200)
    I = sinc_intensity_profile(x, res["a"], 578e-9, D, 1.0)
    plt.plot(x * 1e3, I)
    plt.xlabel("x (mm)")
    plt.ylabel("Intensity (a.u.)")
    plt.title(f"sinc check (a={res['a']*1e6:.0f} um)")
    plt.grid(True, alpha=0.3)
    plt.show()