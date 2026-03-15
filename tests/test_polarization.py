import numpy as np
from optics import polarization

def test_fit_malus_law_runs():
    angles = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    intensities = np.array([91.7, 91.3, 85.7, 75.7, 60.7, 44.5, 29.2, 16.8, 12.8, 0.0])

    result = polarization.fit_malus_law(angles, intensities)

    assert "I0" in result
    assert "offset" in result
    assert "R2" in result
    assert result["I0"] > 0
    assert 0 <= result["R2"] <= 1