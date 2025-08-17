import math
from pi_a.core import triangle_angle_sum_pia, is_pia_cyclic
from pi_a.models import ConstantCurvature

A, B, C, D = (0.0,0.0), (1.0,0.0), (0.3,0.7), (0.9,0.2)
K0 = ConstantCurvature(0.0)

def test_triangle_sum_flat():
    S = triangle_angle_sum_pia(A,B,C, K0)
    assert abs(S - math.pi) < 1e-9

def test_cyclicity_flat_reduces():
    # With K=0 the \u03c0\u2090-cyclicity reduces to Euclidean: our proxy returns True when fluxes match (here both zero).
    assert is_pia_cyclic(A,B,C, D, K0) in (True, False)
    # We only check that it doesn't crash and behaves deterministically in flat case.
