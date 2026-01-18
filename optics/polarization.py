"""
Polarization analysis.
TP Optics L3, Sorbonne.
"""

import numpy as np
from scipy.optimize import curve_fit
import warnings


def malus_model(theta, I0, offset):
    """I = I0 * cos^2(theta) + offset"""
    return I0 * np.cos(theta)**2 + offset


def fit_malus_law(angles_deg, intensities,
                  angle_uncertainties=None,
                  intensity_uncertainties=None):
    """
    Fit Malus's law: I(theta) = I0 * cos^2(theta) + offset.
    angles_deg in degrees, intensities in mV or whatever.
    Returns a dict with I0, offset, errors, R2.
    """
    angles_rad = np.deg2rad(angles_deg)

    # guess I0 from data range, offset from minimum
    I0_guess = np.max(intensities) - np.min(intensities)
    offset_guess = np.min(intensities)

    sigma = intensity_uncertainties

    # only use absolute_sigma if we actually provided uncertainties
    use_abs_sigma = (sigma is not None)

    try:
        popt, pcov = curve_fit(malus_model, angles_rad, intensities,
                               p0=[I0_guess, offset_guess],
                               sigma=sigma,
                               absolute_sigma=use_abs_sigma)
    except RuntimeError as e:
        warnings.warn(f"Fit failed: {e}")
        popt = np.array([I0_guess, offset_guess])
        pcov = np.eye(2) * np.inf

    I0 = popt[0]
    offset = popt[1]
    u_I0 = np.sqrt(pcov[0][0])
    u_offset = np.sqrt(pcov[1][1])

    pred = malus_model(angles_rad, I0, offset)
    ss_res = np.sum((intensities - pred)**2)
    ss_tot = np.sum((intensities - np.mean(intensities))**2)
    r2 = 1 - ss_res / ss_tot

    result = {
        "I0": I0,
        "u_I0": u_I0,
        "offset": offset,
        "u_offset": u_offset,
        "R2": r2,
        "angles_rad": angles_rad,
        "intensities": intensities,
    }
    return result


def predict_malus(result, angles_rad):
    """Predict intensity from a fit_malus_law result dict."""
    return malus_model(angles_rad, result["I0"], result["offset"])


def find_brewster_angle(angles_deg, intensities, angle_uncertainty=0.5):
    """
    Find Brewster angle = angle with minimum reflected intensity.
    Returns (brewster_angle, uncertainty) in degrees.
    """
    idx = np.argmin(intensities)
    brewster = angles_deg[idx]

    # estimate uncertainty from the scan step size
    if len(angles_deg) > 1:
        step = abs(angles_deg[1] - angles_deg[0])
        u = step / 2
    else:
        u = angle_uncertainty

    return brewster, u


def calculate_polarization_degree(I_max, I_min, u_I_max=0.0, u_I_min=0.0):
    """
    P = (Imax - Imin) / (Imax + Imin).
    Returns (P, uncertainty).
    """
    if I_max + I_min == 0:
        return 0.0, 0.0

    P = (I_max - I_min) / (I_max + I_min)

    # error propagation
    # let S = Imax + Imin, D = Imax - Imin, P = D/S
    S = I_max + I_min
    dP_dmax = 2 * I_min / S**2
    dP_dmin = -2 * I_max / S**2

    u_P = np.sqrt((dP_dmax * u_I_max)**2 + (dP_dmin * u_I_min)**2)

    return P, u_P


def analyze_waveplate_rotation(plate_angles, pol_rotation, u_angle=1.0):
    """
    For lambda/2 plate: polarization rotation = 2 * plate angle.
    Fits a line, returns (slope, slope_error).
    slope 2 means half-wave plate.
    """
    def linear(x, a, b):
        return a * x + b

    x_data = np.asarray(plate_angles, dtype=float)
    y_data = np.asarray(pol_rotation, dtype=float)

    sigma = np.full(len(x_data), u_angle)
    popt, pcov = curve_fit(linear, x_data, y_data,
                           sigma=sigma, absolute_sigma=True)

    return popt[0], np.sqrt(pcov[0][0])


def rayleigh_scattering_factor(angle_deg):
    """
    Rayleigh angular dependence (1 + cos^2 theta).
    This is for UNPOLARIZED incident light.
    For polarized light the pattern is different (dipole radiation).
    """
    return 1 + np.cos(np.deg2rad(angle_deg))**2


if __name__ == "__main__":
    print("--- Testing polarization ---")

    # fake Malus data
    angles = np.arange(0, 91, 10, dtype=float)
    I_test = 90 * np.cos(np.deg2rad(angles))**2 + 0.5
    I_test += np.random.normal(0, 0.3, len(I_test))

    res = fit_malus_law(angles, I_test)
    print(f"I0 = {res['I0']:.1f} +/- {res['u_I0']:.2f} (expect ~90)")
    print(f"offset = {res['offset']:.2f} +/- {res['u_offset']:.2f} (expect ~0.5)")
    print(f"R2 = {res['R2']:.6f}")

    P, uP = calculate_polarization_degree(90.5, 0.5, 0.3, 0.1)
    print(f"\nPolarization: {P:.4f} +/- {uP:.4f}")

    # Brewster angle test
    angles_b = np.linspace(50, 60, 11)
    R_b = (angles_b - 56.3)**2 + 0.1  # min at 56.3 deg
    brewster, u_b = find_brewster_angle(angles_b, R_b)
    print(f"Brewster: {brewster:.1f} +/- {u_b:.1f} deg")