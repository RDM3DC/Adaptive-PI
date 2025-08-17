from pi_a.core import polygon_angle_sum_pia
from pi_a.models import ConstantCurvature

square = [(0,0),(1,0),(1,1),(0,1)]

def test_polygon_sum_flat_square():
    S = polygon_angle_sum_pia(square, ConstantCurvature(0.0))
    # Sum of interior angles of a square = 2\u03c0
    import math
    assert abs(S - 2*math.pi) < 1e-9
