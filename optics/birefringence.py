"""
Birefringence analysis.
TP Optics L3, Sorbonne.
"""

import warnings
import numpy as np
from scipy.optimize import curve_fit
from scipy import signal


def phase_difference(dn, e, lam):
    """dphi = 2*pi*dn*e / lambda"""
    return 2 * np.pi * dn * e / lam


def analyze_spectral_fringes(wavelengths, intensities, e, smoothing=True):
    """
    Get birefringence from spectral fringes (cannelures).
    Fringes are equally spaced in 1/lambda, spacing = 1/(dn*e).
    Returns (dn, list of peak wavelengths).
    """
    if smoothing:
        # savgol window — tried a few values, this works for our spectrometer
        wlen = min(51, len(intensities) // 10)
        if wlen % 2 == 0:
            wlen += 1
        from scipy.signal import savgol_filter
        intensities_smooth = savgol_filter(intensities, wlen, 3)
    else:
        intensities_smooth = intensities

    # these peak-finding params were tuned by hand for our TP spectra
    peaks, _ = signal.find_peaks(intensities_smooth,
                                 prominence=0.1 * np.max(intensities),
                                 distance=len(wavelengths) // 20)

    if len(peaks) < 2:
        warnings.warn("Not enough peaks")
        return np.nan, []

    fringe_wl = wavelengths[peaks]

    # spacing in 1/lambda - dn*e
    inv_lam = 1 / fringe_wl
    spacings = np.diff(inv_lam)
    dn = 1 / (np.mean(spacings) * e)

    return dn, fringe_wl.tolist()


def stress_birefringence(sigma, C=None):
    """
    dn = C * sigma (stress-optic law).
    C defaults to 3.5e-12 for PMMA (from Hecht table 8.1 / our TP sheet).
    Change it if you're using glass or polycarbonate.
    """
    if C is None:
        C = 3.5e-12  # PMMA
    return C * sigma


def wave_plate_type(dn, e, lam):
    """
    Check what kind of wave plate this is.
    Returns (retardation_in_waves, type_string).
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


def find_extinction_angles(angles, intensities, threshold_frac=0.1):
    """
    Find angles where intensity is near zero (extinction).
    Returns list of angles.
    """
    threshold = threshold_frac * np.max(intensities)

    minima, _ = signal.find_peaks(-intensities,
                                  height=-threshold,
                                  distance=len(angles) // 10)

    return angles[minima].tolist()


def check_circular_polarization(I_0, I_90, I_45, I_135, tol=0.1):
    """
    If intensity is about the same at all 4 angles -> circular.
    Returns (is_circular, uniformity).
    """
    vals = np.array([I_0, I_90, I_45, I_135])
    mean_I = np.mean(vals)
    if mean_I == 0:
        return False, 0.0

    variation = np.std(vals) / mean_I
    return variation < tol, 1 - variation


def malus_through_birefringent(theta_deg, dn, e, lam):
    """
    Intensity through birefringent plate between crossed polarizers.
    I = sin^2(2*theta) * sin^2(dphi/2)
    """
    theta = np.deg2rad(theta_deg)
    dphi = 2 * np.pi * dn * e / lam
    return np.sin(2 * theta)**2 * np.sin(dphi / 2)**2


if __name__ == "__main__":
    print(" Testing birefringence ")

    # quartz numbers
    dn = 0.009
    e = 1e-3   # 1 mm
    lam = 589e-9  # sodium D

    dphi = phase_difference(dn, e, lam)
    print(f"Phase diff: {dphi/np.pi:.1f} pi rad")

    ret, ptype = wave_plate_type(dn, e, lam)
    print(f"Retardation: {ret:.2f} waves -> {ptype}")

    # simulate fringes and recover dn
    wls = np.linspace(400e-9, 700e-9, 2000)
    ph = 2 * np.pi * dn * e / wls
    I = np.sin(ph)**2

    dn_recovered, peaks = analyze_spectral_fringes(wls, I, e)
    print(f"dn recovered = {dn_recovered:.4f} (expected {dn})")
    print(f"Found {len(peaks)} peaks")