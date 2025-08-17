from pi_a.models import ConstantCurvature, GaussianBump
from pi_a.core import sector_flux
import math

# --- Power-of-a-Point under πₐ -------------------------------------------

def power_of_point_pia(P, A, B, C, D, K):
    """Compute adaptive power of point P wrt circle through A,B,C,D.

    Euclidean version: PA·PB = PC·PD when A,B,C,D concyclic.
    πₐ version: insert flux-based correction factor.
    """
    PA = math.dist(P,A); PB = math.dist(P,B)
    PC = math.dist(P,C); PD = math.dist(P,D)
    base = PA*PB - PC*PD
    # Flux correction: approximate imbalance around P wrt arcs
    Phi1 = sector_flux(P, A, B, K)
    Phi2 = sector_flux(P, C, D, K)
    Xi = 0.5*(Phi1 - Phi2)
    return base * math.exp(Xi)

# --- Simple demo -----------------------------------------------------------
A,B,C,D = (0.0,0.0),(2.0,0.0),(1.0,1.732),(1.0,-1.732)
P = (1.0,0.5)

print("Flat space power:", power_of_point_pia(P,A,B,C,D, ConstantCurvature(0.0)))
print("Curved space power:", power_of_point_pia(P,A,B,C,D, GaussianBump(0.02,1.0,0.5,0.6)))
