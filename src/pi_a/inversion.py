from __future__ import annotations
from typing import Tuple
import math
import numpy as np
from .models import CurvatureField

Point = Tuple[float,float]

class Circle:
    def __init__(self, O: Point, r: float):
        self.O, self.r = O, r


def tangent_sector_flux(P: Point, C: Circle, K: CurvatureField, samples: int = 256) -> float:
    """Approximate curvature flux over the (Euclidean) lens between circle and its tangents from P.

    Construction (flat proxy):
      - If P is outside: draw the two tangent points T1,T2 from P to circle.
      - Region decomposes into two narrow sectors at those tangency points; we approximate
        total flux by sampling along the arc they span plus linear variation toward P.
      - If P is inside, we use a symmetric wedge approximation using diameter through P.

    This is a *first-order* correction: for small |K| and moderate distances it stabilizes
    power comparisons under mild curvature.
    """
    (px,py), (ox,oy), r = P, C.O, C.r
    dx, dy = px-ox, py-oy
    d2 = dx*dx + dy*dy
    if d2 == 0:
        return 0.0
    d = math.sqrt(d2)
    if d < r:  # interior point: construct pseudo-tangent via orthogonal diameter
        # Choose orthogonal direction to (dx,dy)
        nx, ny = -dy/d, dx/d
        # Two boundary points along diameter through direction (nx,ny)
        T1 = (ox + r*nx, oy + r*ny)
        T2 = (ox - r*nx, oy - r*ny)
    else:
        # Outside: classic tangent points
        # angle between OP and OT: theta = arccos(r/d)
        theta = math.acos(min(1.0, max(-1.0, r/d)))
        base_ang = math.atan2(dy, dx)
        a1 = base_ang + theta
        a2 = base_ang - theta
        T1 = (ox + r*math.cos(a1), oy + r*math.sin(a1))
        T2 = (ox + r*math.cos(a2), oy + r*math.sin(a2))

    # Sample along arc from T1 to T2 (choose shorter arc through interior)
    aT1 = math.atan2(T1[1]-oy, T1[0]-ox)
    aT2 = math.atan2(T2[1]-oy, T2[0]-ox)
    # unwrap shortest
    da = aT2 - aT1
    if da <= -math.pi:
        da += 2*math.pi
    elif da > math.pi:
        da -= 2*math.pi
    # discretize arc & radial interpolation toward P
    m = max(8, samples//8)
    arc_pts = []
    for i in range(m+1):
        a = aT1 + da * (i / m)
        x = ox + r*math.cos(a)
        y = oy + r*math.sin(a)
        arc_pts.append((x,y))
    # For each arc sample, integrate curvature along segment to P (one-point mid rule)
    flux = 0.0
    for (x,y) in arc_pts:
        mx = 0.5*(x+px)
        my = 0.5*(y+py)
        # weight ~ area of small triangle (O,x,P) via 0.5*| (x-O) x (P-O) | / m as proxy
        # cross product magnitude between (x-ox,y-oy) and (px-ox,py-oy)
        cx = (x-ox)*(py-oy) - (y-oy)*(px-ox)
        tri_area = 0.5*abs(cx) / m
        flux += K(mx,my) * tri_area
    return float(flux)


def power_of_point_pia(P: Point, C: Circle, K: CurvatureField) -> float:
    """πₐ power of a point: Euclidean power plus curvature lens flux correction.

    Power(P;C) = |OP|^2 - r^2  (Euclidean)  + Φ_lens  (first-order curvature term)
    """
    (px,py), (ox,oy), r = P, C.O, C.r
    euclid = (px-ox)**2 + (py-oy)**2 - r**2
    Phi = tangent_sector_flux(P, C, K)
    return euclid + Phi


def invert_point_pia(P: Point, C: Circle, K: CurvatureField) -> Point:
    (x,y), (ox,oy), r2 = P, C.O, C.r**2
    dx, dy = x-ox, y-oy
    d2 = dx*dx + dy*dy
    if d2 == 0:
        raise ValueError("Cannot invert center")
    s = r2 / d2
    return (ox + s*dx, oy + s*dy)


def power_decomposition_pia(P: Point, C: Circle, K: CurvatureField):
    """Return (euclidean_power, curvature_flux) components."""
    (x,y), (ox,oy), r = P, C.O, C.r
    euclid = (x-ox)**2 + (y-oy)**2 - r**2
    flux = tangent_sector_flux(P, C, K)
    return euclid, flux
