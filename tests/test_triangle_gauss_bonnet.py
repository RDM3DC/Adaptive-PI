import math
from pi_a.core import triangle_angle_sum_pia
from pi_a.models import ConstantCurvature
from pi_a.geometry import tri_area


def test_gauss_bonnet_triangle_sum():
    A, B, C = (0.0, 0.0), (1.0, 0.0), (0.0, 1.0)
    K = ConstantCurvature(0.05)
    S = triangle_angle_sum_pia(A, B, C, K)
    expected = math.pi + 0.05 * tri_area(A, B, C)
    assert abs(S - expected) < 1e-2
