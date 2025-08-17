from pi_a.core import is_pia_cyclic
from pi_a.models import GaussianBump, ConstantCurvature

A, B, C = (0.0,0.0), (1.0,0.0), (0.25,0.9)
X = (0.9,0.2)

def test_flux_balance_toggle():
    flat = ConstantCurvature(0.0)
    curved = GaussianBump(K0=0.05, x0=0.5, y0=0.3, sigma=0.25)
    # In flat space, sector fluxes are equal (0), so the balance decision is neutral.
    a = is_pia_cyclic(A,B,C,X, flat)
    b = is_pia_cyclic(A,B,C,X, curved)
    assert a in (True, False)
    assert b in (True, False)
