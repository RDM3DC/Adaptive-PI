from pi_a.inversion import Circle, power_of_point_pia, power_decomposition_pia, invert_point_pia
from pi_a.models import ConstantCurvature, GaussianBump

C = Circle((0.0, 0.0), 1.0)
P = (2.0, 0.0)

# Flat case
base = power_of_point_pia(P, C, ConstantCurvature(0.0))
print("Euclidean power (should be 3):", base)
print("Inversion of P across C:", invert_point_pia(P, C, ConstantCurvature(0.0)))

# Curved case (first-order placeholder still near base)
_, corr = power_decomposition_pia(P, C, GaussianBump(K0=0.02, x0=0.3, y0=0.2, sigma=0.5))
print("πₐ correction (placeholder, expect ~0 for now):", corr)
