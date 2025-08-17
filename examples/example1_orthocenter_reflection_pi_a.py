import math
from pi_a.core import triangle_angle_sum_pia, is_pia_cyclic
from pi_a.models import ConstantCurvature, GaussianBump

A, B, C = (0.0, 0.0), (1.2, 0.0), (0.25, 0.9)

# Simplified orthocenter reflection demo: we only check \u03c0\u2090-cyclicity of (A,B,C,H')
# Here we fabricate an H' near the circumcircle to illustrate flux\u2011balanced behavior.
Hprime = (0.95, 0.18)

print("Flat case (K=0):")
K0 = ConstantCurvature(0.0)
print("  cyclic? ", is_pia_cyclic(A,B,C,Hprime, K0))

print("Curved case (Gaussian bump near center):")
Kb = GaussianBump(K0=+0.02, x0=0.5, y0=0.35, sigma=0.3)
print("  cyclic? ", is_pia_cyclic(A,B,C,Hprime, Kb))
