from pi_a.core import is_pia_cyclic
from pi_a.models import ConstantCurvature, GaussianBump

# Triangle with orthocenter H
A, B, C = (0.0, 0.0), (2.0, 0.0), (0.5, 1.5)
# Hardcode an approximate reflection of orthocenter across BC for demo
Hprime = (0.9, 0.2)

print("Euclidean test (K=0):", is_pia_cyclic(A,B,C, Hprime, ConstantCurvature(0.0)))
print("πₐ test (Gaussian bump):", is_pia_cyclic(A,B,C, Hprime, GaussianBump(K0=0.02, x0=0.7, y0=0.6, sigma=0.5)))
