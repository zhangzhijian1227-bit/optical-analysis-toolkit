# Quick Start

## Installation
```bash
git clone https://github.com/zhangzhijian1227-bit/optical-analysis-toolkit.git
cd optical-analysis-toolkit
pip install -r requirements.txt
```

Check that it works:
```bash
python -c "from optics import polarization; print('ok')"
```

## Basic examples

### Malus's law
```python
from optics import polarization
import numpy as np

angles = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
intensities = np.array([91.7, 91.3, 85.7, 75.7, 60.7, 44.5, 29.2, 16.8, 12.8, 0.0])

res = polarization.fit_malus_law(angles, intensities)
print(res['I0'], res['offset'], res['R2'])

# predict at some angle
I_45 = polarization.predict_malus(res, np.deg2rad(np.array([45])))
```

### Single slit diffraction
```python
from optics import diffraction
import numpy as np

positions = np.array([-1.35, -0.93, -0.42, 0.45, 0.89, 1.36]) * 1e-3
orders = np.array([-3, -2, -1, 1, 2, 3])

res = diffraction.analyze_single_slit(positions, orders, D=0.30, lam=578e-9)
print(f"a = {res['a']*1e6:.1f} +/- {res['a_err']*1e6:.1f} um")
```

### Birefringence
```python
from optics import birefringence
import numpy as np

dn = 0.009
e = 2e-3
lam = 589e-9

dphi = birefringence.phase_difference(dn, e, lam)
print(f"Phase difference: {dphi/np.pi:.2f} pi rad")

ret, ptype = birefringence.wave_plate_type(dn, e, lam)
print(f"Plate type: {ptype}")
```

### With uncertainties
```python
res = polarization.fit_malus_law(
    angles, intensities,
    intensity_uncertainties=np.full_like(intensities, 0.1)
)
```

## Notebooks

The `examples/` folder has full analysis notebooks on real lab data — `malus_law_demo.ipynb` is probably the best starting point.

```bash
cd examples
jupyter notebook malus_law_demo.ipynb
```

All functions have docstrings (`help(polarization.fit_malus_law)`). If something breaks, open an issue or email zhangzhijian1227@gmail.com.