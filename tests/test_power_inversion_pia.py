from pi_a.inversion import Circle, power_of_point_pia, power_decomposition_pia, invert_point_pia
from pi_a.models import ConstantCurvature, GaussianBump
import math

C = Circle((0.0, 0.0), 1.0)
P = (2.0, 0.0)

def test_power_flat():
    assert abs(power_of_point_pia(P, C, ConstantCurvature(0.0)) - 3.0) < 1e-9

def test_inversion_flat_roundtrip():
    Q = invert_point_pia(P, C, ConstantCurvature(0.0))  # should be (0.5, 0.0)
    assert abs(Q[0] - 0.5) < 1e-9 and abs(Q[1]) < 1e-12

def test_power_correction_placeholder():
    base, corr = power_decomposition_pia(P, C, GaussianBump(K0=0.02, x0=0.3, y0=0.2, sigma=0.5))
    assert abs(base - 3.0) < 1e-9
    assert abs(corr) < 1e-6  # placeholder is ~0
