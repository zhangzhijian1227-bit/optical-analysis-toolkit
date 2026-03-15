import numpy as np
from optics import birefringence

def test_phase_difference():
    dn = 0.01
    thickness = 1e-3
    wl = 500e-9
    expected = 2 * np.pi * dn * thickness / wl
    result = birefringence.phase_difference(dn, thickness, wl)
    assert np.isclose(result, expected)