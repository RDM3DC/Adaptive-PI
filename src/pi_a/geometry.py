from __future__ import annotations
import math
from typing import Tuple

Point = Tuple[float, float]

def tri_area(A: Point, B: Point, C: Point) -> float:
    return 0.5 * ((B[0]-A[0])*(C[1]-A[1]) - (B[1]-A[1])*(C[0]-A[0]))

def angle_at(P: Point, Q: Point, R: Point) -> float:
    """Euclidean interior angle \u2220PQR in radians."""
    u = (P[0]-Q[0], P[1]-Q[1])
    v = (R[0]-Q[0], R[1]-Q[1])
    nu = math.hypot(*u)
    nv = math.hypot(*v)
    if nu == 0 or nv == 0:
        return 0.0
    dot = (u[0]*v[0] + u[1]*v[1])/(nu*nv)
    dot = max(-1.0, min(1.0, dot))
    return math.acos(dot)
