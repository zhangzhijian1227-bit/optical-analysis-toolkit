import pytest
from optics import birefringence

def test_phase_difference_zero_wavelength():
    with pytest.raises(ZeroDivisionError):
        birefringence.phase_difference(0.01, 1e-3, 0)