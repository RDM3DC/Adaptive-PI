from __future__ import annotations
from typing import Iterable, Tuple, List
import math
from .models import CurvatureField

Point = Tuple[float,float]

def geodesic_circle(O: Point, r: float, K: CurvatureField, n: int = 128) -> List[Point]:
    """Sample a 'geodesic circle' using Euclidean proxy; for small curvature treat as Euclidean circle.
    Future refinement: adjust radial distance based on average curvature along rays.
    """
    pts = []
    for i in range(n):
        th = 2*math.pi * i / n
        pts.append((O[0] + r*math.cos(th), O[1] + r*math.sin(th)))
    return pts
