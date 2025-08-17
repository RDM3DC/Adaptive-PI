from __future__ import annotations
from typing import Tuple
from .models import CurvatureField


class Circle:
    def __init__(self, O: Tuple[float,float], r: float):
        self.O, self.r = O, r


def power_of_point_pia(P, C: Circle, K: CurvatureField) -> float:
    (x,y), (ox,oy), r = P, C.O, C.r
    euclid = (x-ox)**2 + (y-oy)**2 - r**2
    return euclid  # TODO: add flux correction


def invert_point_pia(P, C: Circle, K: CurvatureField):
    (x,y), (ox,oy), r2 = P, C.O, C.r**2
    dx, dy = x-ox, y-oy
    d2 = dx*dx + dy*dy
    if d2 == 0:
        raise ValueError("Cannot invert center")
    s = r2 / d2
    return (ox + s*dx, oy + s*dy)
