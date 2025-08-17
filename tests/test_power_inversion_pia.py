from pi_a.inversion import Circle, power_of_point_pia, invert_point_pia
from pi_a.models import ConstantCurvature

C = Circle((0,0), 1.0)
K0 = ConstantCurvature(0.0)

def test_power_center():
    P = (0.0,0.0)
    val = power_of_point_pia(P,C,K0)
    assert abs(val - (-1.0)) < 1e-9

def test_invert_outside():
    P = (2.0,0.0)
    Q = invert_point_pia(P,C,K0)
    assert abs(Q[0]-0.5) < 1e-9
