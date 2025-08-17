from pi_a.inversion import Circle, power_of_point_pia, invert_point_pia, tangent_sector_flux
from pi_a.models import ConstantCurvature, GaussianBump
import math

C = Circle((0,0), 1.0)
K0 = ConstantCurvature(0.0)

def test_power_center_flat():
    P = (0.0,0.0)
    val = power_of_point_pia(P,C,K0)
    # Center power = -r^2
    assert abs(val - (-1.0)) < 1e-9

def test_invert_outside_flat():
    P = (2.0,0.0)
    Q = invert_point_pia(P,C,K0)
    assert abs(Q[0]-0.5) < 1e-9

def test_double_inversion_is_identity():
    C2 = Circle((0.3, -0.4), 1.7)
    P = (2.2, 0.9)
    K_flat = ConstantCurvature(0.0)
    Q = invert_point_pia(P, C2, K_flat)
    P2 = invert_point_pia(Q, C2, K_flat)
    assert math.hypot(P2[0]-P[0], P2[1]-P[1]) < 1e-9

def test_curvature_flux_nonzero():
    bump = GaussianBump(K0=0.05, x0=0.2, y0=-0.1, sigma=0.4)
    P = (1.5, 0.6)
    flat = power_of_point_pia(P, C, K0)
    curved = power_of_point_pia(P, C, bump)
    # flux may increase or decrease; just ensure difference captured
    assert abs(curved - flat) > 0.0

def test_tangent_sector_flux_symmetry_on_axes():
    # For P on x-axis, lens symmetry; flux should be finite and stable
    P = (2.0, 0.0)
    val = tangent_sector_flux(P, C, K0)
    assert abs(val) < 1e-12  # flat curvature gives zero
