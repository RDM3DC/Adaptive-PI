import math
from pi_a.core import triangle_angle_sum_pia
from pi_a.models import ConstantCurvature, GaussianBump
from pi_a.geometry import angle_at

# This example just prints Euclidean vs \u03c0\u2090 angle sums as a prelude to Ceva-trig forms.
A, B, C = (0.0, 0.0), (1.0, 0.0), (0.2, 0.8)

K0 = ConstantCurvature(0.0)
Kb = GaussianBump(K0=0.03, x0=0.4, y0=0.3, sigma=0.35)

S_flat = triangle_angle_sum_pia(A,B,C, K0)
S_curv = triangle_angle_sum_pia(A,B,C, Kb)

print("Angle sum (flat):   ", S_flat)
print("Angle sum (curved): ", S_curv, "  (\u0394=", S_curv - S_flat, ")")

# NOTE: Trig\u2011Ceva with explicit flux\u2011correction factor is documented in SPEC; a
# numerical demo would require explicit cevian split points and sector flux terms.
