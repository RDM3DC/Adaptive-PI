from pi_a.inversion import Circle, power_of_point_pia, invert_point_pia
from pi_a.models import ConstantCurvature

# Simple IMO-style demonstration: Given circle (O,r) and external point P, compute power and inverse point.

C = Circle((0.0,0.0), 3.0)
K = ConstantCurvature(0.0)
P = (7.0, 1.0)

powP = power_of_point_pia(P, C, K)
print(f"Power of P wrt C: {powP:.6f}")

invP = invert_point_pia(P, C, K)
print(f"Inversion of P: {invP}")
