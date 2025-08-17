from pi_a.ceva import trig_ceva_pia, concurrent_pia
from pi_a.models import ConstantCurvature, GaussianBump

# Triangle and a concurrent triple (medians -> centroid)
A, B, C = (0.0, 0.0), (2.0, 0.0), (0.5, 1.8)
D = ((B[0]+C[0])/2.0, (B[1]+C[1])/2.0)  # on BC
E = ((C[0]+A[0])/2.0, (C[1]+A[1])/2.0)  # on CA
F = ((A[0]+B[0])/2.0, (A[1]+B[1])/2.0)  # on AB

K_flat = ConstantCurvature(0.0)
K_bump = GaussianBump(K0=0.02, x0=0.7, y0=0.6, sigma=0.5)

def test_trig_ceva_flat_is_one():
    val = trig_ceva_pia(A,B,C, D,E,F, K_flat, use_flux=False)
    assert abs(val - 1.0) < 1e-9

def test_concurrent_pia_flat():
    assert concurrent_pia(A,B,C, D,E,F, K_flat, tol=1e-6)

def test_concurrent_pia_curved_tolerant():
    # Under curvature, the corrected value should remain ~1 for symmetric medians.
    assert concurrent_pia(A,B,C, D,E,F, K_bump, tol=5e-2)

def test_nonconcurrent_breaks():
    # Move F off the midpoint to break concurrency
    F_bad = (F[0] + 0.3, F[1])
    val = trig_ceva_pia(A,B,C, D,E,F_bad, K_flat, use_flux=False)
    assert abs(val - 1.0) > 1e-3
