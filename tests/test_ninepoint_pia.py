from pi_a.ninepoint import ninepoint_curve_pia
from pi_a.models import ConstantCurvature, GaussianBump

A,B,C = (0.0,0.0),(2.0,0.0),(0.5,1.7)

def test_ninepoint_flat_exists():
    O,R,pts = ninepoint_curve_pia(A,B,C, ConstantCurvature(0.0))
    assert all(isinstance(p, tuple) for p in pts)

def test_ninepoint_curved_runs():
    Kb = GaussianBump(K0=0.03, x0=0.6, y0=0.5, sigma=0.5)
    O,R,pts = ninepoint_curve_pia(A,B,C, Kb)
    assert O and R >= 0
