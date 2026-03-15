# Optical Analysis Toolkit

A Python toolkit for analyzing optics lab data (polarization, diffraction, birefringence), developed during L3 Physics at Sorbonne Université.

## Why does it exist?

Basically it's a set of Python functions to process data from optics labs. Instead of redoing everything by hand in Excel every time, I put all the fits, uncertainty calculations and plots together in one package.

The data used in the examples comes from actual lab sessions:
- Polarization (Malus' law, waveplates, Brewster's angle)
- Diffraction (single/double slit, gratings)
- Birefringence (channeled spectra, photoelasticity)

## What it actually does

**`polarization.py`** - Malus' law fits ($I = I_0 \cos^2\theta$), waveplate analysis ($\lambda$/2, $\lambda$/4), Brewster's angle, degree of polarization. Returns fitted parameters, uncertainties and goodness-of-fit $R^2$.

**`diffraction.py`** - Extracts slit widths using $sinc^2$ profile fitting for single slit, fringe spacing for double slits, and grating characterization.

**`birefringence.py`** - Calculates $\Delta n = \Delta/e$, analyzes channeled spectra (spectral fringes), and handles photoelasticity data.

**`visualization.py`** - Generates report-ready plots with error bars, residuals, polar plots.

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
print(f"Slit width: {result.slit_width*1e6:.1f} +/- {result.u_slit_width*1e6:.1f} µm")
```

### Birefringence
```python
from optics import birefringence

delta_n = birefringence.calculate_birefringence(
    optical_path_difference=1.2e-6,
    thickness=0.5e-3
)
```


## Notebooks

The notebooks in `examples/` show how to use the toolkit on real lab data. The most complete one is `malus_law.ipynb` which does the full Malus' law analysis from start to finish.

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
