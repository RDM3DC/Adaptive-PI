from __future__ import annotations
import math
from typing import Tuple
from .geometry import angle_at
from .core import sector_flux
from .models import CurvatureField

Point = Tuple[float,float]

def opposite_angle_sum_pia(A:Point,B:Point,C:Point,D:Point, K:CurvatureField, center:Point) -> float:
    sum_e = angle_at(A,B,C) + angle_at(A,D,C)
    Phi = 0.5*( sector_flux(center, B, C, K) - sector_flux(center, D, A, K) )
    return sum_e + Phi


def is_pia_cyclic_quad(A:Point,B:Point,C:Point,D:Point, K:CurvatureField, center:Point, tol:float=1e-3) -> bool:
    S = opposite_angle_sum_pia(A,B,C,D, K, center)
    return abs(S - math.pi) < tol
