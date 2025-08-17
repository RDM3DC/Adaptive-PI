from __future__ import annotations
import math
from typing import Tuple
from .geometry import angle_at
from .core import sector_flux
from .models import CurvatureField

Point = Tuple[float,float]

def trig_menelaus_pia(A: Point, B: Point, C: Point,
                                                                                        P: Point, Q: Point, R: Point,
                                                                                        K: CurvatureField,
                                                                                        use_flux: bool = True) -> float:
                """Flux‑corrected trig Menelaus identity using vertex-based formulation.

                To avoid degenerate π angles at transversal points lying on lines, we use
                equivalent vertex expression (see standard references):
                        (sin ∠ABP / sin ∠PBC) * (sin ∠BCQ / sin ∠QCA) * (sin ∠CAR / sin ∠RAB) = 1
                Provided P∈AB, Q∈BC, R∈CA.
                """
                ABP = angle_at(A,B,P); PBC = angle_at(P,B,C)
                BCQ = angle_at(B,C,Q); QCA = angle_at(Q,C,A)
                CAR = angle_at(C,A,R); RAB = angle_at(R,A,B)
                num = math.sin(ABP) * math.sin(BCQ) * math.sin(CAR)
                den = math.sin(PBC) * math.sin(QCA) * math.sin(RAB)
                base = num/den if den != 0 else float('inf')
                if not use_flux:
                                return base
                # Flux corrections: difference between numerator and denominator sector pairs at each transversal point
                Phi_P = sector_flux(P, A, B, K) - sector_flux(P, C, B, K)
                Phi_Q = sector_flux(Q, B, C, K) - sector_flux(Q, A, C, K)
                Phi_R = sector_flux(R, C, A, K) - sector_flux(R, B, A, K)
                Xi = 0.5*(Phi_P + Phi_Q + Phi_R)
                return base * math.exp(Xi)

def collinear_pia(A: Point, B: Point, C: Point,
                  P: Point, Q: Point, R: Point,
                  K: CurvatureField, tol: float = 1e-3) -> bool:
    val = trig_menelaus_pia(A,B,C, P,Q,R, K, use_flux=True)
    return math.isfinite(val) and abs(val - 1.0) < tol
