import math
from pi_a.ceva import trig_ceva_pia, concurrent_pia
from pi_a.models import ConstantCurvature

A,B,C = (0.0,0.0), (1.0,0.0), (0.3,0.7)
D,E,F = (0.5,0.0), (0.65,0.35), (0.15,0.35)

K0 = ConstantCurvature(0.0)

def test_trig_ceva_flat():
    val = trig_ceva_pia(A,B,C, D,E,F, K0)
    assert abs(val-1.0) < 1e-7

def test_concurrent_flat():
    assert concurrent_pia(A,B,C, D,E,F, K0)
