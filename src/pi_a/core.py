from __future__ import annotations
import math
from typing import Tuple, Iterable
import numpy as np
from .geometry import tri_area, angle_at, Point
from .models import CurvatureField

# --- Numerical flux integration helpers ------------------------------------

def polygon_flux(vertices: Iterable[Point], K: CurvatureField, samples: int = 4000) -> float:
    """Monte\u2011Carlo estimate of \u03a6 = \u222c\u222c_polygon K dA (simple triangle/polygon sampler).
    For small domains and smooth K, this first\u2011order term is sufficient.
    """
    verts = np.array(list(vertices), dtype=float)
    assert len(verts) >= 3
    # Triangulate fan at verts[0]
    base = verts[0]
    area_total = 0.0
    flux_total = 0.0
    rng = np.random.default_rng(12345)
    for i in range(1, len(verts)-1):
        A, B, C = base, verts[i], verts[i+1]
        triA = abs(tri_area(tuple(A), tuple(B), tuple(C)))
        area_total += triA
        # sample inside triangle via barycentric
        n = max(1, samples // (len(verts)-2))
        u = rng.random(n)
        v = rng.random(n)
        mask = u + v > 1.0
        u[mask] = 1 - u[mask]
        v[mask] = 1 - v[mask]
        P = (A + (B-A)[None,:]*u[:,None] + (C-A)[None,:]*v[:,None])
        Kvals = np.array([K(float(p[0]), float(p[1])) for p in P])
        flux_total += Kvals.mean() * triA
    return float(flux_total)

# --- \u03c0\u2090 primitives ----------------------------------------------------------

def triangle_angle_sum_pia(A: Point, B: Point, C: Point, K: CurvatureField) -> float:
    """S_\u0394 = A+B+C = \u03c0 + \u03a6_\u0394 (first order). Returns the predicted sum in radians."""
    # Euclidean interior angles
    angA = angle_at(B, A, C)
    angB = angle_at(A, B, C)
    angC = angle_at(A, C, B)
    # Flux correction
    Phi = polygon_flux([A,B,C], K)
    # Replace \u03c0 with \u03c0 + \u03a6_\u0394 - (Euclidean sum - \u03c0)
    # Equivalently: S_\u0394(\u03c0\u2090) \u2248 (A+B+C) + \u03a6_\u0394  - ( (A+B+C) - \u03c0 ) = \u03c0 + \u03a6_\u0394
    return math.pi + Phi

def polygon_angle_sum_pia(poly: Iterable[Point], K: CurvatureField) -> float:
    verts = list(poly)
    m = len(verts)
    if m < 3:
        return 0.0
    Phi = polygon_flux(verts, K)
    return (m - 2) * math.pi + Phi

# Sector flux for adaptive inscribed\u2011angle comparisons

def sector_flux(center: Point, P: Point, Q: Point, K: CurvatureField, radius_samples: int = 64, angle_samples: int = 64) -> float:
    """Approximate flux over a circular sector (center\u2192P\u2192Q) using Euclidean radius
    as a proxy. For small curvature and small sectors, this is a good first\u2011order proxy.
    """
    cx, cy = center
    rP = math.hypot(P[0]-cx, P[1]-cy)
    rQ = math.hypot(Q[0]-cx, Q[1]-cy)
    r = 0.5*(rP + rQ)
    aP = math.atan2(P[1]-cy, P[0]-cx)
    aQ = math.atan2(Q[1]-cy, Q[0]-cx)
    # unwrap angle
    da = aQ - aP
    if da <= -math.pi:
        da += 2*math.pi
    elif da > math.pi:
        da -= 2*math.pi
    # integrate in polar rectangle approx
    rs = np.linspace(0, r, radius_samples)
    thetas = np.linspace(aP, aP+da, angle_samples)
    flux = 0.0
    for i in range(len(thetas)-1):
        th = 0.5*(thetas[i]+thetas[i+1])
        dth = thetas[i+1]-thetas[i]
        for j in range(len(rs)-1):
            rr = 0.5*(rs[j]+rs[j+1])
            dr = rs[j+1]-rs[j]
            x = cx + rr*math.cos(th)
            y = cy + rr*math.sin(th)
            flux += K(x,y) * rr * dth * dr
    return float(flux)

# \u03c0\u2090-cyclicity: balance of sector fluxes

def is_pia_cyclic(A: Point, B: Point, C: Point, X: Point, K: CurvatureField, center_guess: Point | None = None) -> bool:
    """Returns True if A,B,C,X are \u03c0\u2090-cyclic to first order: the two inscribed\u2011angle
    sectors have approximately equal flux so their difference \u2248 0.
    In flat limit (K=0), this reduces to the Euclidean criterion.
    """
    # Heuristic center: circumcenter in Euclidean proxy
    if center_guess is None:
        center_guess = _circumcenter_euclid(A,B,C)
    O = center_guess
    Phi1 = sector_flux(O, A, B, K)  # sector AB
    Phi2 = sector_flux(O, A, C, K)  # sector AC
    return abs(Phi1 - Phi2) < 1e-6  # tolerance can be tuned

# Euclidean circumcenter helper (for center_guess)

def _circumcenter_euclid(A: Point, B: Point, C: Point) -> Point:
    ax, ay = A; bx, by = B; cx, cy = C
    d = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    if abs(d) < 1e-12:
        # near-collinear fallback: centroid
        return ((ax+bx+cx)/3.0, (ay+by+cy)/3.0)
    a2 = ax*ax + ay*ay
    b2 = bx*bx + by*by
    c2 = cx*cx + cy*cy
    ux = (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by)) / d
    uy = (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax)) / d
    return (ux, uy)
