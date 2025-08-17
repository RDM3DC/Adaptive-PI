from __future__ import annotations
from typing import Tuple, List
import math
from .models import CurvatureField
from .circumcurve import circumcurve_pia_balanced

Point = Tuple[float,float]

def _line(P:Point,Q:Point):
    (x1,y1),(x2,y2)=P,Q; return (y1-y2, x2-x1, x1*y2-x2*y1)

def _intersect(L1,L2):
    A1,B1,C1=L1; A2,B2,C2=L2
    D=A1*B2-A2*B1
    return ((B1*C2-B2*C1)/D, (C1*A2-C2*A1)/D)

def _perp_through(L, P):
    A,B,C=L; x0,y0=P; return (B,-A, -(B*x0 - A*y0))

def orthocenter(A:Point,B:Point,C:Point) -> Point:
    BC=_line(B,C); CA=_line(C,A)
    altA=_perp_through(BC,A); altB=_perp_through(CA,B)
    return _intersect(altA,altB)

def nine_points(A:Point,B:Point,C:Point) -> List[Point]:
    H = orthocenter(A,B,C)
    mid = lambda P,Q: ((P[0]+Q[0])/2.0,(P[1]+Q[1])/2.0)
    D,E,F = mid(B,C), mid(C,A), mid(A,B)
    def foot(P, L):
        A,B,C=L; x0,y0=P; t=(A*x0+B*y0+C)/(A*A+B*B); return (x0 - A*t, y0 - B*t)
    BC=_line(B,C); CA=_line(C,A); AB=_line(A,B)
    Ha,Hb,Hc = foot(A,BC), foot(B,CA), foot(C,AB)
    Da,Db,Dc = mid(A,H), mid(B,H), mid(C,H)
    return [D,E,F, Ha,Hb,Hc, Da,Db,Dc]

def ninepoint_curve_pia(A:Point,B:Point,C:Point, K:CurvatureField):
    pts = nine_points(A,B,C)
    O,R,_ = circumcurve_pia_balanced(pts[0],pts[1],pts[2], K)
    return O,R,pts
