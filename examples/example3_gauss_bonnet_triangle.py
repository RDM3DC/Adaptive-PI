"""Demonstrate Gauss--Bonnet correction for triangle angle sums."""
import math
from pi_a.core import triangle_angle_sum_pia
from pi_a.models import ConstantCurvature
from pi_a.geometry import tri_area

A, B, C = (0.0, 0.0), (1.0, 0.0), (0.0, 1.0)
K = ConstantCurvature(0.05)  # constant positive curvature

S = triangle_angle_sum_pia(A, B, C, K)
expected = math.pi + 0.05 * tri_area(A, B, C)

print("Angle sum with curvature:", S)
print("Expected (\u03c0 + flux):", expected)
print("Deviation from flat \u03c0:", S - math.pi)
