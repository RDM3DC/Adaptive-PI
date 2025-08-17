from pi_a.projective import cross_ratio, cross_ratio_pia

A,B = (0.0,0.0),(2.0,0.0)
C,D = (0.5,0.2),(1.4,-0.15)

def test_cross_ratio_affine_invariance():
    A2,B2 = (1.0,1.0),(3.0,1.0)
    C2,D2 = (1.5,1.2),(2.4,0.85)
    assert abs(cross_ratio(A,B,C,D) - cross_ratio(A2,B2,C2,D2)) < 1e-12

def test_cross_ratio_pia_first_order():
    assert abs(cross_ratio_pia(A,B,C,D) - cross_ratio(A,B,C,D)) < 1e-12
