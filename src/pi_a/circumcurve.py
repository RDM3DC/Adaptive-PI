from __future__ import annotations
from typing import Tuple
import math
from .geometry import tri_area

Point = Tuple[float,float]

def euclid_circumcenter(A: Point, B: Point, C: Point) -> Point:
    ax, ay = A; bx, by = B; cx, cy = C
    d = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    if abs(d) < 1e-12:
        return ((ax+bx+cx)/3.0, (ay+by+cy)/3.0)
    a2 = ax*ax + ay*ay; b2 = bx*bx + by*by; c2 = cx*cx + cy*cy
    ux = (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by)) / d
    uy = (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax)) / d
    return (ux, uy)

def euclid_circumradius(A: Point, B: Point, C: Point) -> float:
    ax, ay = A; bx, by = B; cx, cy = C
    a = math.hypot(bx-cx, by-cy)
    b = math.hypot(ax-cx, ay-cy)
    c = math.hypot(ax-bx, ay-by)
    area = abs(tri_area(A,B,C))
    if area == 0:
        return 0.0
    return (a*b*c) / (4*area)
