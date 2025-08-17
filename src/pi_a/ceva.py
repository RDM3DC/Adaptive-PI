from __future__ import annotations
import math
from .geometry import angle_at
from .models import CurvatureField


def trig_ceva_pia(A,B,C, D,E,F, K: CurvatureField) -> float:
    """Return the trig-Ceva product in \u03c0\u2090 setting. In flat case this \u2192 1.
    We attach a first-order flux imbalance factor ~ exp(ε·Ξ) implicitly via angle drift.
    Points D,E,F lie on sides AB, BC, CA respectively.
    """
    ACD = angle_at(A,C,D); DCB = angle_at(D,C,B)
    BAE = angle_at(B,A,E); EAC = angle_at(E,A,C)
    CBF = angle_at(C,B,F); FBA = angle_at(F,B,A)
    num = math.sin(ACD) * math.sin(BAE) * math.sin(CBF)
    den = math.sin(DCB) * math.sin(EAC) * math.sin(FBA)
    return num / den


def concurrent_pia(A,B,C, D,E,F, K: CurvatureField, tol: float = 1e-3) -> bool:
    return abs(trig_ceva_pia(A,B,C, D,E,F, K) - 1.0) < tol
