"""
Diffraction analysis.
TP Optics L3, Sorbonne.
"""

import numpy as np
from scipy.optimize import curve_fit

def linear_func(x, a, b):
    return a * x + b


def analyze_single_slit(x_positions, orders, D, lam=578e-9):
    """
    Calculate slit width 'a'.
    Formula: x = (lambda * D / a) * p
    """
    popt, pcov = curve_fit(linear_func, orders, x_positions)

    slope = popt[0]  # = lam * D / a
    intercept = popt[1]
    slope_err = np.sqrt(pcov[0][0])

    # a = lam * D / slope
    a = (lam * D) / slope
    a_err = a * (slope_err / slope)

    residuals = x_positions - linear_func(orders, slope, intercept)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((x_positions - np.mean(x_positions))**2)
    r2 = 1 - (ss_res / ss_tot)

    results = {
        "a": a,
        "a_err": a_err,
        "R2": r2,
        "slope": slope,
        "intercept": intercept,
    }
    return results


def analyze_double_slit(x_positions, orders, D, lam=578e-9):
    """
    Calculate slit separation 'b'.
    Same idea: x = (lambda * D / b) * n
    """
    popt, pcov = curve_fit(linear_func, orders, x_positions)

    slope = popt[0]
    slope_err = np.sqrt(pcov[0][0])

    b = (lam * D) / slope
    b_err = b * (slope_err / slope)

    return b, b_err


def get_wavelength_grating(angles_deg, orders, d):
    """
    Find wavelength from grating with known d.
    d * sin(theta) = n * lambda
    => sin(theta) = (lambda/d) * n
    """
    theta = np.deg2rad(angles_deg)
    y_vals = np.sin(theta)

    def grating_model(n, slope):
        return slope * n

    popt, pcov = curve_fit(grating_model, orders, y_vals)

    slope = popt[0]  # = lambda / d
    slope_err = np.sqrt(pcov[0][0])

    wl = slope * d
    wl_err = slope_err * d

    return wl, wl_err


def get_grating_constant(sin_theta, wavelengths):
    """
    Find grating constant d from known spectral lines (1st order).
    sin(theta) = (1/d) * lambda
    Used this for the mercury lamp in TP5.
    """
    def model(lam, A):
        return A * lam

    popt, pcov = curve_fit(model, wavelengths, sin_theta)

    A = popt[0]  # = 1/d
    A_err = np.sqrt(pcov[0][0])

    d = 1 / A
    d_err = A_err / A**2

    return d, d_err


def sinc_intensity_profile(x, a, lam, f, I0):
    """
    Theoretical single slit intensity.
    I = I0 * (sin(beta)/beta)^2
    beta = (k * a * x) / (2 * f)
    """
    k = 2 * np.pi / lam
    beta = (k * a * x) / (2 * f)

    # add tiny number to avoid division by zero
    # make it works
    beta = beta + 1e-10

    intensity = I0 * (np.sin(beta) / beta)**2
    return intensity


def fringe_visibility(I_max, I_min):
    """V = (Imax - Imin) / (Imax + Imin)."""
    if I_max + I_min == 0:
        return 0.0
    return (I_max - I_min) / (I_max + I_min)


if __name__ == "__main__":
    print(" Testing diffraction ")

    # fake data
    orders = np.array([-3, -2, -1, 1, 2, 3])
    D = 0.3  # 30 cm
    pos = np.array([-0.0013, -0.00085, -0.00042, 0.00044, 0.00088, 0.00132])

    res = analyze_single_slit(pos, orders, D)
    print(f"a = {res['a']*1e6:.1f} um")
    print(f"error = {res['a_err']*1e6:.1f} um")
    print(f"R2 = {res['R2']:.6f}")

    # check sinc shape
    import matplotlib.pyplot as plt
    x = np.linspace(-0.005, 0.005, 200)
    I = sinc_intensity_profile(x, res["a"], 578e-9, D, 1.0)
    plt.plot(x * 1e3, I)
    plt.xlabel("Position x (mm)")
    plt.ylabel("Intensity (a.u.)")
    plt.title(f"Sinc check (a={res['a']*1e6:.0f}um)")
    plt.grid(True, alpha=0.3)
    plt.show()