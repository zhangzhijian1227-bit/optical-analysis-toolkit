"""
Polarization analysis.
TP Optics L3, Sorbonne.
"""

import numpy as np
from scipy.optimize import curve_fit


def malus_model(theta, I0, offset):
    return I0 * np.cos(theta)**2 + offset

def linear_func(x, a, b):
    return a * x + b


def fit_malus_law(angles_deg, intensities, sigma=None):
    """
    Fit Malus's law: I(theta) = I0 * cos^2(theta) + offset.
    Returns dict with I0, offset, errors, R2.
    """
    angles_rad = np.deg2rad(angles_deg)

    I0_guess = np.max(intensities) - np.min(intensities)
    offset_guess = np.min(intensities)

    try:
        popt, pcov = curve_fit(malus_model, angles_rad, intensities,
                               p0=[I0_guess, offset_guess],
                               sigma=sigma, absolute_sigma=(sigma is not None))
    except RuntimeError:
        print("fit_malus: failed to converge")
        return None

    I0, offset = popt
    u_I0 = np.sqrt(pcov[0][0])
    u_offset = np.sqrt(pcov[1][1])

    pred = malus_model(angles_rad, I0, offset)
    ss_res = np.sum((intensities - pred)**2)
    ss_tot = np.sum((intensities - np.mean(intensities))**2)
    r2 = 1 - ss_res / ss_tot

    return {
        "I0": I0, "u_I0": u_I0,
        "offset": offset, "u_offset": u_offset,
        "R2": r2,
        "angles_rad": angles_rad,
        "intensities": intensities,
    }


def find_brewster_angle(angles, intensities):
    """
    Brewster angle = angle of minimum reflected intensity.
    Uncertainty estimated from angular step size.
    """
    idx = np.argmin(intensities)
    brewster = angles[idx]
    step = abs(angles[1] - angles[0]) if len(angles) > 1 else 0.5
    return brewster, step / 2


def polarization_degree(I_max, I_min, u_Imax=0.0, u_Imin=0.0):
    """P = (Imax - Imin) / (Imax + Imin), with error propagation."""
    if I_max + I_min == 0:
        return 0.0, 0.0

    P = (I_max - I_min) / (I_max + I_min)
    S = I_max + I_min
    # error propagation
    da = 2 * I_min / S**2
    db = -2 * I_max / S**2
    u_P = np.sqrt((da * u_Imax)**2 + (db * u_Imin)**2)

    return P, u_P


def analyze_waveplate_rotation(phi_plate, phi_pol, u_angle=1.0):
    """
    lambda/2 plate: polarization rotation = 2 * plate angle.
    Linear fit, slope ~ 2 for a half-wave plate.
    """
    x = np.asarray(phi_plate, dtype=float)
    y = np.asarray(phi_pol, dtype=float)
    sigma = np.full(len(x), u_angle)

    popt, pcov = curve_fit(linear_func, x, y, sigma=sigma, absolute_sigma=True)
    return popt[0], np.sqrt(pcov[0][0])


def rayleigh_scattering_factor(angle_deg):
    # Rayleigh scattering, unpolarized incident light
    return 1 + np.cos(np.deg2rad(angle_deg))**2


if __name__ == "__main__":
    print("--- testing polarization ---")

    angles = np.arange(0, 91, 10, dtype=float)
    I_test = 90 * np.cos(np.deg2rad(angles))**2 + 0.5
    I_test += np.random.normal(0, 0.3, len(I_test))

    res = fit_malus_law(angles, I_test)
    if res:
        print(f"I0 = {res['I0']:.1f} +/- {res['u_I0']:.2f}  (expect ~90)")
        print(f"offset = {res['offset']:.2f} +/- {res['u_offset']:.2f}  (expect ~0.5)")
        print(f"R2 = {res['R2']:.6f}")

    P, uP = polarization_degree(90.5, 0.5, 0.3, 0.1)
    print(f"\nPolarization degree: {P:.4f} +/- {uP:.4f}")

    angles_b = np.linspace(50, 60, 11)
    R_b = (angles_b - 56.3)**2 + 0.1
    brewster, u_b = find_brewster_angle(angles_b, R_b)
    print(f"Brewster: {brewster:.1f} +/- {u_b:.1f} deg")