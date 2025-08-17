from pi_a.ceva import trig_ceva_pia, concurrent_pia
from pi_a.models import ConstantCurvature, GaussianBump

A, B, C = (0.0, 0.0), (2.0, 0.0), (0.5, 1.8)
D = ((B[0]+C[0])/2.0, (B[1]+C[1])/2.0)
E = ((C[0]+A[0])/2.0, (C[1]+A[1])/2.0)
F = ((A[0]+B[0])/2.0, (A[1]+B[1])/2.0)

print("Euclidean trig-Ceva product:", trig_ceva_pia(A,B,C, D,E,F, ConstantCurvature(0.0), use_flux=False))
print("πₐ concurrency:", concurrent_pia(A,B,C, D,E,F, GaussianBump(K0=0.02, x0=0.7, y0=0.6, sigma=0.5)))
