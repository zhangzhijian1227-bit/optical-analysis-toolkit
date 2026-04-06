"""
Plotting helpers for optics experiments.
TP Optics L3, Sorbonne.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# stub copies to avoid relative import issues when running standalone
# in the real project: from .polarization import malus_model
def malus_model(theta, I0, offset):
    return I0 * np.cos(theta)**2 + offset

def sinc_intensity_profile(x, a, lam, f, I0):
    X = (a * x) / (lam * f)
    return I0 * np.sinc(X)**2


def set_plot_style():
    """Call this at the start to get consistent TP-style plots."""
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


def plot_malus_law(angles, I, fit, u_angles=None, u_I=None,
                   title="Malus's law", save_path=None):
    """
    Plot Malus's law data + fit.
    fit is the dict from fit_malus_law.
    """
    I0 = fit['I0']
    offset = fit['offset']

    fig, ax = plt.subplots()

    ax.errorbar(angles, I, xerr=u_angles, yerr=u_I,
                fmt='o', capsize=4, label='data', alpha=0.8)

    theta_fine = np.linspace(np.deg2rad(angles.min()), np.deg2rad(angles.max()), 200)
    ax.plot(np.rad2deg(theta_fine), malus_model(theta_fine, I0, offset),
            'r-', lw=2, label=f'fit  $I_0$={I0:.1f}')

    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Intensity (mV)')
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def plot_diffraction(x_pos, I, a, lam, f, title="Diffraction pattern", save_path=None):
    """
    Measured diffraction pattern with sinc^2 envelope.
    x_pos in meters, a = slit width, lam = wavelength, f = focal length.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(x_pos * 1e3, I, 'o', ms=5, alpha=0.8, label='measured')

    span = x_pos.max() - x_pos.min()
    x_th = np.linspace(x_pos.min() - 0.1*span, x_pos.max() + 0.1*span, 500)
    ax.plot(x_th * 1e3, sinc_intensity_profile(x_th, a, lam, f, np.max(I)),
            'r-', lw=2, label='sinc$^2$', alpha=0.7)

    ax.set_xlabel('x (mm)')
    ax.set_ylabel('Intensity (a.u.)')
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def plot_spectrum(wl, I, peaks=None, title="Spectrum", save_path=None):
    """
    Plot spectrum with optional peak markers.
    wl in nm, peaks is a list of wavelengths in nm.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(wl, I, '-', lw=1.5)

    if peaks is not None and len(peaks) > 0:
        peaks = np.asarray(peaks)
        if peaks.mean() < 1e-3:   # passed in meters by mistake
            peaks = peaks * 1e9
        ax.plot(peaks, np.interp(peaks, wl, I), 'ro', ms=8, label='peaks')
        ax.legend()

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Intensity (counts)')
    ax.set_title(title)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    set_plot_style()

    # Malus
    angles = np.arange(0, 361, 20, dtype=float)
    I_noisy = 80 * np.cos(np.deg2rad(angles))**2 + 5 + np.random.normal(0, 2, len(angles))
    plot_malus_law(angles, I_noisy, {'I0': 80, 'offset': 5},
                   u_I=np.full(len(angles), 2))

    # diffraction
    x_pos = np.linspace(-0.01, 0.01, 100)
    I_diff = sinc_intensity_profile(x_pos, 100e-6, 632.8e-9, 1.0, 100)
    I_diff += np.random.normal(0, 1.5, len(x_pos))
    I_diff = np.clip(I_diff, 0, None)
    plot_diffraction(x_pos, I_diff, 100e-6, 632.8e-9, 1.0,
                     title="Single slit (a=100 um)")

    # spectrum
    wl = np.linspace(400, 700, 500)
    ph = 2 * np.pi * 0.01 * 0.001 / (wl * 1e-9)
    I_spec = np.sin(ph)**2 * 1000 + np.random.normal(0, 50, len(wl))
    peak_idx, _ = find_peaks(I_spec, prominence=200)
    plot_spectrum(wl, I_spec, peaks=wl[peak_idx])
 
    plt.show()