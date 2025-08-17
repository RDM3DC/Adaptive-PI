from __future__ import annotations
from typing import Tuple
import math
from .geometry import tri_area
from .core import sector_flux
from .models import CurvatureField

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
    b = math.hypot(ax-cx, ay-cx)
    c = math.hypot(ax-bx, ay-by)
    area = abs(tri_area(A,B,C))
    if area == 0:
        return 0.0
    return (a*b*c) / (4*area)

def circumcurve_pia(A: Point, B: Point, C: Point, K: CurvatureField):
    """First-order πₐ circumcurve: returns Euclidean circumcenter & radius.
    Placeholder for curvature-adjusted center (future balancing step)."""
    O = euclid_circumcenter(A,B,C)
    R = math.hypot(A[0]-O[0], A[1]-O[1])
    return O, R

def circumcurve_pia_balanced(A: Point, B: Point, C: Point, K: CurvatureField, iters: int = 12, step: float = 0.15):
    """Attempt to balance sector flux differences at the three vertices by nudging center.
    Returns (O,R,residual) where residual is sum of absolute flux imbalances.
    Current implementation: small gradient-free heuristic; in flat limit reduces to Euclidean.
    """
    Ox,Oy = euclid_circumcenter(A,B,C)
    def residual(Ox,Oy):
        PhiA = sector_flux((Ox,Oy), B, C, K) - sector_flux((Ox,Oy), C, A, K)
        PhiB = sector_flux((Ox,Oy), C, A, K) - sector_flux((Ox,Oy), A, B, K)
        PhiC = sector_flux((Ox,Oy), A, B, K) - sector_flux((Ox,Oy), B, C, K)
        return abs(PhiA) + abs(PhiB) + abs(PhiC)
    base = residual(Ox,Oy)
    # Simple coordinate descent perturbations
    h = 1e-3
    for _ in range(iters):
        base = residual(Ox,Oy)
        rx = residual(Ox+h, Oy)
        if rx < base:
            Ox += step*h
            continue
        ry = residual(Ox, Oy+h)
        if ry < base:
            Oy += step*h
            continue
        step *= 0.5
        if step < 1e-4:
            break
    R = math.hypot(A[0]-Ox, A[1]-Oy)
    return (Ox,Oy), R, residual(Ox,Oy)
