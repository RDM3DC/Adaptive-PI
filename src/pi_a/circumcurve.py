from __future__ import annotations
from typing import Tuple
from .geodesics import geodesic_circle
from .models import CurvatureField
import math

Point = Tuple[float, float]

def circumcurve_pia(A: Point, B: Point, C: Point, K: CurvatureField) -> tuple[Point, float]:
    """Return an adaptive 'circumcurve' center and radius for points A,B,C.

    v0: Euclidean circumcenter/radius (flat-limit exact). Later versions will
    iteratively nudge the center to balance sector-flux differences.
    """
    (ax,ay),(bx,by),(cx,cy) = A,B,C
    d = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    if abs(d) < 1e-12:
        # Degenerate: return centroid and NaN radius
        return ((ax+bx+cx)/3.0, (ay+by+cy)/3.0), float('nan')
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux = (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by))/d
    uy = (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax))/d
    R = math.hypot(ax-ux, ay-uy)
    return (ux,uy), R
