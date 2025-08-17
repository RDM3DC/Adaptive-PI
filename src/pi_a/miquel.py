from __future__ import annotations
from typing import Tuple, List
import math
from .core import sector_flux
from .models import CurvatureField

Point = Tuple[float,float]

def _circumcenter(P: Point, Q: Point, R: Point) -> Tuple[Point, float]:
    (ax,ay),(bx,by),(cx,cy) = P,Q,R
    d = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    if abs(d) < 1e-12:
        return ((ax+bx+cx)/3.0,(ay+by+cy)/3.0), float('nan')
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux = (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by))/d
    uy = (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax))/d
    R = math.hypot(ax-ux, ay-uy)
    return (ux,uy), R

def _circle_intersections(O1: Point, R1: float, O2: Point, R2: float) -> List[Point]:
    (x1,y1),(x2,y2) = O1,O2
    dx,dy = x2-x1, y2-y1
    d = math.hypot(dx,dy)
    if d == 0: return []
    a = (R1*R1 - R2*R2 + d*d) / (2*d)
    h2 = R1*R1 - a*a
    if h2 < 0: return []
    xm = x1 + a*dx/d
    ym = y1 + a*dy/d
    rx = -dy/d
    ry = dx/d
    h = math.sqrt(max(0.0,h2))
    return [(xm + h*rx, ym + h*ry), (xm - h*rx, ym - h*ry)]

def miquel_point_euclid(A: Point, B: Point, C: Point, D: Point) -> Point:
    O1,R1 = _circumcenter(A,B,C)
    O2,R2 = _circumcenter(B,C,D)
    cand = _circle_intersections(O1,R1, O2,R2)
    if not cand:
        return ((A[0]+B[0]+C[0]+D[0])/4.0, (A[1]+B[1]+C[1]+D[1])/4.0)
    O3,R3 = _circumcenter(C,D,A)
    def resid(P):
        return abs(math.hypot(P[0]-O3[0], P[1]-O3[1]) - R3)
    return min(cand, key=resid)

def miquel_point_pia(A: Point, B: Point, C: Point, D: Point, K: CurvatureField,
                     iters: int = 20, step: float = 0.1, fd_eps: float = 1e-3) -> Tuple[Point, float]:
    Mx,My = miquel_point_euclid(A,B,C,D)
    def triple_residual(M: Point) -> float:
        def res(triple):
            O,_ = _circumcenter(*triple)
            A1,B1,C1 = triple
            return abs(sector_flux(O, A1, B1, K) - sector_flux(O, A1, C1, K))
        return res((A,B,M)) + res((B,C,M)) + res((C,D,M)) + res((D,A,M))
    for _ in range(iters):
        S = triple_residual((Mx,My))
        Sx = triple_residual((Mx+fd_eps,My))
        Sy = triple_residual((Mx,My+fd_eps))
        Mx -= step * (Sx - S)/fd_eps
        My -= step * (Sy - S)/fd_eps
    return (Mx,My), triple_residual((Mx,My))
