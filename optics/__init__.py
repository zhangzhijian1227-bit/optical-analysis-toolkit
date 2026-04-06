"""
Optical Analysis Toolkit
L3 Physics, Sorbonne Université.
"""

__version__ = "1.0.0"
__author__ = "ZHANG ZHIJIAN"

from . import polarization
from . import diffraction
from . import birefringence
from . import visualization

from .polarization import fit_malus_law, find_brewster_angle
from .diffraction import analyze_single_slit, analyze_double_slit
from .birefringence import phase_difference, analyze_spectral_fringes
from .visualization import plot_malus_law, plot_diffraction, plot_spectrum