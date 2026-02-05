# Optical Analysis Toolkit

Python toolkit for analyzing optics lab data (polarization, diffraction, birefringence). Made for L3 Physics at Sorbonne Université.

## What is this?

Basically it's a set of Python functions to process data from optics labs. Instead of redoing everything by hand in Excel every time, I put all the fits, uncertainty calculations and plots together in one package.

The data used in the examples comes from actual lab sessions:
- Polarization (Malus' law, waveplates, Brewster's angle)
- Diffraction (single/double slit, gratings)
- Birefringence (channeled spectra, photoelasticity)

## Structure
```
optics/
├── polarization.py    # Malus' law, waveplates, Brewster
├── diffraction.py     # Single/double slit, gratings
├── birefringence.py   # Birefringence, channeled spectra
└── visualization.py   # Plotting functions

examples/
├── malus_law.ipynb
├── diffraction_analysis.ipynb
└── birefringence_analysis.ipynb
```

## Installation

You need Python 3.7+ and the usual dependencies:
```bash
git clone https://github.com/zhangzhijian1227-bit/optical-analysis-toolkit.git
cd optical-analysis-toolkit
pip install -r requirements.txt
```

Or as a package:
```bash
pip install -e .
```

## Usage

### Polarization - Malus' Law
```python
from optics import polarization, visualization
import numpy as np

angles = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
intensities = np.array([91.7, 91.3, 85.7, 75.7, 60.7, 44.5, 29.2, 16.8, 12.8, 0.0])

result = polarization.fit_malus_law(angles, intensities)
print(result)

fig, ax = visualization.plot_malus_law(angles, intensities, result)
```

### Diffraction
```python
from optics import diffraction

# Determine slit width from the diffraction profile
result = diffraction.fit_single_slit(positions, intensities, wavelength=633e-9)
print(f"Slit width: {result.slit_width*1e6:.1f} ± {result.u_slit_width*1e6:.1f} µm")
```

### Birefringence
```python
from optics import birefringence

delta_n = birefringence.calculate_birefringence(
    optical_path_difference=1.2e-6,
    thickness=0.5e-3
)
```

## What it actually does

**`polarization.py`** - Fits Malus' law ($I = I_0 \cos^2\theta$), waveplate analysis (λ/2, λ/4), Brewster's angle determination, degree of polarization. Fits return parameters with uncertainties and R².

**`diffraction.py`** - sinc² profile fitting for single slit, fringe spacing for double slit, grating characterization. You can also compute theoretical profiles to compare with.

**`birefringence.py`** - Birefringence calculation Δn = Δ/e, channeled spectra analysis (spectral fringes), waveplate classification, photoelasticity.

**`visualization.py`** - Plots with error bars, residuals, polar plots for polarization, spectral plots. The style is clean enough to put in a lab report.

## Notebooks

The notebooks in `examples/` show how to use the toolkit on real lab data. The most complete one is `malus_law_demo.ipynb` which does the full Malus' law analysis from start to finish.

## Dependencies

- numpy >= 1.19
- scipy >= 1.5
- matplotlib >= 3.3
- jupyter (for the notebooks)

## Author

ZHANG Zhijian
L3 Physics, Sorbonne Université

## License

MIT - see LICENSE file.
