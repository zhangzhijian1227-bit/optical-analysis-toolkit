"""
Plotting helpers for optics experiments.
TP Optics L3, Sorbonne.
"""

import numpy as np
import matplotlib.pyplot as plt

from .polarization import malus_model
from .diffraction import sinc_intensity_profile


def set_plot_style():
    """
    Set some nicer defaults for matplotlib.
    Call this yourself if you want the TP-style plots.
    """
    plt.rcParams.update({
        'font.size': 11,
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'errorbar.capsize': 3,
        'figure.figsize': (9, 6),
        'axes.grid': True,
        'grid.alpha': 0.3,
    })


def plot_malus_law(angles_deg, intensities, fit_result,
                   angle_errors=None, intensity_errors=None,
                   title="Malus's law", save_path=None):
    """
    Plot Malus's law data + fit + residuals.
    fit_result is the dict from fit_malus_law.
    Returns (fig, axes).
    """
    I0 = fit_result['I0']
    offset = fit_result['offset']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]},
                                    sharex=True)

    # data
    ax1.errorbar(angles_deg, intensities,
                 xerr=angle_errors, yerr=intensity_errors,
                 fmt='o', capsize=4, label='data', alpha=0.8)

    # fit curve
    angles_fine_rad = np.linspace(np.deg2rad(min(angles_deg)),
                                   np.deg2rad(max(angles_deg)), 200)
    I_fit = malus_model(angles_fine_rad, I0, offset)
    ax1.plot(np.rad2deg(angles_fine_rad), I_fit, 'r-', lw=2, label='fit')

    ax1.set_ylabel('Intensity')
    ax1.set_title(title)
    ax1.legend()
    ax1.grid(alpha=0.3)

    # residuals
    pred = malus_model(np.deg2rad(angles_deg), I0, offset)
    residuals = intensities - pred

    ax2.errorbar(angles_deg, residuals, yerr=intensity_errors,
                 fmt='o', capsize=3, alpha=0.8)
    ax2.axhline(0, color='red', ls='--', lw=1)
    ax2.set_xlabel('Angle (deg)')
    ax2.set_ylabel('Residual')
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig, (ax1, ax2)


def plot_diffraction_pattern(positions, intensities, a, lam, f,
                             title="Diffraction pattern", save_path=None):
    """
    Plot measured diffraction pattern with theoretical sinc^2 envelope.
    positions in meters, a = slit width, lam = wavelength, f = focal length.
    Returns (fig, ax).
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    # data
    ax.plot(positions * 1e3, intensities, 'o',
            label='measured', ms=5, alpha=0.8)

    # theory — extend range a bit so the curve doesn't get cut off
    span = positions.max() - positions.min()
    x_th = np.linspace(positions.min() - 0.1 * span,
                        positions.max() + 0.1 * span, 500)
    I_th = sinc_intensity_profile(x_th, a, lam, f, np.max(intensities))
    ax.plot(x_th * 1e3, I_th, 'r-', lw=2, label='sinc^2 theory', alpha=0.7)

    ax.set_xlabel('Position (mm)')
    ax.set_ylabel('Intensity')
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig, ax


def plot_spectral_fringes(wavelengths_nm, intensities,
                          peak_wavelengths=None,
                          title="Spectral fringes", save_path=None):
    """
    Plot spectrum with optional peak markers.
    peak_wavelengths can be in meters or nm — we auto-detect.
    Returns (fig, ax).
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(wavelengths_nm, intensities, '-', lw=1.5)

    if peak_wavelengths is not None and len(peak_wavelengths) > 0:
        peaks_nm = np.array(peak_wavelengths)
        # if values are tiny they're probably in meters, convert to nm
        if np.mean(peaks_nm) < 1e-3:
            peaks_nm = peaks_nm * 1e9

        # interpolate to get the intensity at each peak wavelength
        peak_I = np.interp(peaks_nm, wavelengths_nm, intensities)
        ax.plot(peaks_nm, peak_I, 'ro', ms=8, label='detected peaks')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Intensity')
    ax.set_title(title)
    ax.grid(alpha=0.3)
    if peak_wavelengths is not None:
        ax.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig, ax