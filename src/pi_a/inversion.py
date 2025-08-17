from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Tuple
from .models import CurvatureField

Point = Tuple[float, float]

@dataclass
class Circle:
    O: Point
    r: float

# --- Euclidean baselines ---------------------------------------------------

def power_of_point_euclid(P: Point, C: Circle) -> float:
    (x,y), (ox,oy) = P, C.O
    return (x-ox)*(x-ox) + (y-oy)*(y-oy) - C.r*C.r

# --- πₐ helpers ------------------------------------------------------------

def power_decomposition_pia(P: Point, C: Circle, K: CurvatureField) -> tuple[float, float]:
    """Return (base_euclid_power, first_order_correction).

    For now, we expose the correction as 0.0 (placeholder) while keeping a clean
    hook for future lens/sector flux. This preserves flat-limit exactness and
    lets downstream code budget curvature effects explicitly.
    """
    base = power_of_point_euclid(P, C)
    correction = 0.0  # TODO: lens/sector flux across circle chords
    return base, correction


def power_of_point_pia(P: Point, C: Circle, K: CurvatureField) -> float:
    base, corr = power_decomposition_pia(P, C, K)
    return base + corr


def invert_point_pia(P: Point, C: Circle, K: CurvatureField) -> Point:
    """Invert point P in circle C. Flat-limit exact; πₐ corrections are higher order.

    Returns the Euclidean inversion, which is the correct zeroth/first-order
    mapping when curvature effects are modeled separately via power corrections.
    """
    (x,y), (ox,oy) = P, C.O
    dx, dy = x-ox, y-oy
    d2 = dx*dx + dy*dy
    if d2 == 0:
        raise ValueError("Cannot invert the center")
    s = (C.r*C.r) / d2
    return (ox + s*dx, oy + s*dy)
