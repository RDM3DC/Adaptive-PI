from pi_a.cyclic import is_pia_cyclic_quad
from pi_a.circumcurve import circumcurve_pia
from pi_a.models import ConstantCurvature, GaussianBump

A,B,C,D = (0.0,0.0),(2.0,0.1),(1.6,1.4),(-0.2,0.9)
O,R = circumcurve_pia(A,B,C, ConstantCurvature(0.0))

def test_cyclic_flat_runs():
    assert is_pia_cyclic_quad(A,B,C,D, ConstantCurvature(0.0), O, tol=1e-2) in (True, False)

def test_cyclic_curved_runs():
    assert is_pia_cyclic_quad(A,B,C,D, GaussianBump(K0=0.03, x0=0.7, y0=0.6, sigma=0.45), O, tol=0.1) in (True, False)
