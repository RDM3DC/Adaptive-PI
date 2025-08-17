from __future__ import annotations
from typing import Tuple

Point = Tuple[float,float]

def _project_scalar(A:Point,B:Point,X:Point) -> float:
    ax,ay=A; bx,by=B; xx,yy=X
    ux,uy = bx-ax, by-ay
    vx,vy = xx-ax, yy-ay
    denom = ux*ux+uy*uy
    t = (ux*vx+uy*vy)/denom
    return t

def cross_ratio(A:Point,B:Point,C:Point,D:Point) -> float:
    tC = _project_scalar(A,B,C)
    tD = _project_scalar(A,B,D)
    return (tC*(1-tD)) / (tD*(1-tC))

def cross_ratio_pia(A:Point,B:Point,C:Point,D:Point) -> float:
    return cross_ratio(A,B,C,D)
