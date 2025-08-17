from pi_a.miquel import miquel_point_euclid, miquel_point_pia
from pi_a.models import ConstantCurvature, GaussianBump

A,B,C,D = (0.0,0.0),(2.0,0.2),(0.5,1.6),(-0.4,0.9)

def test_miquel_flat_matches():
    M = miquel_point_euclid(A,B,C,D)
    M2,_ = miquel_point_pia(A,B,C,D, ConstantCurvature(0.0))
    assert abs(M[0]-M2[0]) < 1e-5 and abs(M[1]-M2[1]) < 1e-5

def test_miquel_curved_residual():
    Kb = GaussianBump(K0=0.03, x0=0.7, y0=0.5, sigma=0.45)
    _, S = miquel_point_pia(A,B,C,D, Kb)
    assert S >= 0.0
