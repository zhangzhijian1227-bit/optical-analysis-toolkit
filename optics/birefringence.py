"""
Birefringence analysis.
TP Optics L3, Sorbonne.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter


def phase_difference(dn, e, lam):
    """dphi = 2*pi*dn*e / lambda"""
    return 2 * np.pi * dn * e / lam


def analyze_spectral_fringes(wl, I, e):
    """
    Calculate birefringence dn from spectral fringes.
    """
    # Adjusted prominence and distance manually by looking at the plotted spectra
    peaks, _ = find_peaks(I, prominence=200, distance=50)

    if len(peaks) < 2:
        print("analyze_spectral_fringes: not enough peaks found")
        return np.nan, peaks

    fringe_wl = wl[peaks]
    spacings = np.diff(1 / fringe_wl)

    # dn = 1 / (Delta(1/lambda) * e)
    dn = 1 / (np.mean(spacings) * e)

    return dn, peaks


def stress_birefringence(sigma, C=None):
    """
    dn = C * sigma (stress-optic law).
    C defaults to 3.5e-12 for PMMA (Hecht table 8.1 / TP sheet).
    Change if using glass or polycarbonate.
    """
    if C is None:
        C = 3.5e-12  # PMMA
    return C * sigma


def wave_plate_type(dn, e, lam):
    """
    Check what kind of wave plate this is.
    Returns (retardation in waves, type string).
    """
    ret = dn * e / lam
    frac = ret % 1.0

    # 0.05 wave tolerance, good enough for lab plates
    if abs(frac - 0.25) < 0.05:
        ptype = "lambda/4"
    elif abs(frac - 0.5) < 0.05:
        ptype = "lambda/2"
    elif abs(frac) < 0.05 or abs(frac - 1.0) < 0.05:
        ptype = "full wave"
    else:
        ptype = f"{ret:.2f} waves (not standard)"

    return ret, ptype


def find_extinction_angles(angles, I, threshold_frac=0.1):
    """Find angles where intensity is near zero (extinction)."""
    threshold = threshold_frac * np.max(I)
    minima, _ = find_peaks(-I, height=-threshold, distance=len(angles) // 10)
    return angles[minima]


def check_circular_polarization(I_vals, tol=0.1):
    """
    If intensity is about the same at all angles -> circular.
    Pass I measured at [0, 45, 90, 135] deg.
    Returns (is_circular, uniformity).
    """
    I_vals = np.asarray(I_vals)
    mean_I = np.mean(I_vals)
    if mean_I == 0:
        return False, 0.0
    variation = np.std(I_vals) / mean_I
    return variation < tol, 1 - variation


def malus_through_birefringent(theta_deg, dn, e, lam):
    """
    I = sin^2(2*theta) * sin^2(dphi/2), between crossed polarizers.
    """
    theta = np.deg2rad(theta_deg)
    dphi = 2 * np.pi * dn * e / lam
    return np.sin(2 * theta)**2 * np.sin(dphi / 2)**2


if __name__ == "__main__":
    print("--- testing birefringence ---")

    # quartz numbers
    dn = 0.009
    e = 1e-3    # 1 mm
    lam = 589e-9  # sodium D

    dphi = phase_difference(dn, e, lam)
    print(f"phase diff: {dphi/np.pi:.1f} pi rad")

    ret, ptype = wave_plate_type(dn, e, lam)
    print(f"retardation: {ret:.2f} waves -> {ptype}")

    # simulate fringes and recover dn
    wls = np.linspace(400e-9, 700e-9, 2000)
    I = 1000 * np.sin(2 * np.pi * dn * e / wls)**2

    dn_recovered, peaks = analyze_spectral_fringes(wls, I, e)
    print(f"dn recovered = {dn_recovered:.4f}  (expected {dn})")
    print(f"found {len(peaks)} peaks")