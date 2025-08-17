from __future__ import annotations
import math
from typing import Iterable, Tuple, List
from .models import CurvatureField

Point = Tuple[float, float]

def geodesic_circle(center: Point, r: float, K: CurvatureField, n: int = 256) -> List[Point]:
    """Return n sample points of a 'geodesic circle' around center with radius r.

    First-order proxy: we return Euclidean samples (flat-limit exact). Future
    work can warp samples by accumulated curvature along radial geodesics.
    """
    cx, cy = center
    pts = []
    for i in range(n):
        th = 2*math.pi * i / n
        pts.append((cx + r*math.cos(th), cy + r*math.sin(th)))
    return pts
