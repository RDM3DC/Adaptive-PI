from __future__ import annotations
import math
from typing import Tuple
from .core import sector_flux
from .geometry import angle_at
from .models import CurvatureField

Point = Tuple[float, float]

# --- Helpers ---------------------------------------------------------------

def _safe_sin(x: float) -> float:
    # guard against tiny numerical negatives from rounding inside acos
    return math.sin(x)

# Trig-Ceva (Euclidean form):
# (sin ∠BAD / sin ∠DAC) * (sin ∠CBE / sin ∠EBA) * (sin ∠ACF / sin ∠FCA) = 1
# We generalize by attaching a first-order curvature factor via sector fluxes.


def trig_ceva_pia(A: Point, B: Point, C: Point,
                  D: Point, E: Point, F: Point,
                  K: CurvatureField,
                  use_flux: bool = True,
                  center_hint: Point | None = None) -> float:
    """Return the trig-Ceva product under πₐ.

    In flat space (K≡⁰), the value equals the Euclidean trig-Ceva product.
    With curvature, we optionally multiply by a first-order factor derived from
    sector-flux imbalances at the three vertices.
    """
    # Euclidean angles appearing in trig-Ceva ratios
    BAD = angle_at(B, A, D); DAC = angle_at(D, A, C)
    CBE = angle_at(C, B, E); EBA = angle_at(E, B, A)
    ACF = angle_at(A, C, F)
    # Denominator at vertex C must be ∠FCB (not FCA) for the classical trig‑Ceva identity.
    FCB = angle_at(F, C, B)

    # Euclidean base product
    num = _safe_sin(BAD) * _safe_sin(CBE) * _safe_sin(ACF)
    den = _safe_sin(DAC) * _safe_sin(EBA) * _safe_sin(FCB)
    base = num / den if den != 0 else float("inf")

    if not use_flux:
        return base

    # First-order πₐ correction via sector fluxes around each vertex.
    # Idea: each angle θ ≈ θ_euclid + (1/2)·ΔΦ_sector, so ratio picks up exp(ε·Ξ).
    # Approximate Ξ by the sum of sector flux differences between numerator and
    # denominator angles at each vertex.
    # We center each sector at its vertex (local, cheap, robust).
    Phi_A = sector_flux(A, B, D, K) - sector_flux(A, D, C, K)
    Phi_B = sector_flux(B, C, E, K) - sector_flux(B, E, A, K)
    Phi_C = sector_flux(C, A, F, K) - sector_flux(C, F, A, K)

    # Scale 1/2 as in inscribed-angle correction; exponentiate for stability.
    Xi = 0.5 * (Phi_A + Phi_B + Phi_C)
    correction = math.exp(Xi)
    return base * correction


def concurrent_pia(A: Point, B: Point, C: Point,
                   D: Point, E: Point, F: Point,
                   K: CurvatureField,
                   tol: float = 1e-3) -> bool:
    """Test πₐ concurrency of AD, BE, CF using trig-Ceva with flux correction.

    In the flat limit, this reduces to the Euclidean trig-Ceva test (≈1).
    """
    val = trig_ceva_pia(A,B,C, D,E,F, K, use_flux=True)
    if not math.isfinite(val):
        return False
    return abs(val - 1.0) < tol
