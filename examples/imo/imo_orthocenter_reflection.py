from pi_a.core import is_pia_cyclic
from pi_a.models import ConstantCurvature, GaussianBump
import math

# --- Basic helpers ---------------------------------------------------------

def line_through(P, Q):
    (x1,y1),(x2,y2)=P,Q
    A = y1 - y2
    B = x2 - x1
    C = x1*y2 - x2*y1
    return (A,B,C)

def intersection(L1, L2):
    A1,B1,C1=L1; A2,B2,C2=L2
    D = A1*B2 - A2*B1
    if abs(D) < 1e-12:
        raise ValueError("Lines nearly parallel")
    x = (B1*C2 - B2*C1)/D
    y = (C1*A2 - C2*A1)/D
    return (x,y)

def perpendicular_through(L, P):
    A,B,C = L
    # Perpendicular to Ax+By+C=0 has direction (A,B) swapped -> (B,-A)
    # Line with normal (B,-A) through P has eq: B x - A y + C' = 0
    x0,y0 = P
    C2 = -(B*x0 - A*y0)
    return (B, -A, C2)

def reflect_point_across_line(P, L):
    # Reflect P across line Ax+By+C=0
    A,B,C = L
    x0,y0 = P
    d = (A*x0 + B*y0 + C) / (A*A + B*B)
    xr = x0 - 2*A*d
    yr = y0 - 2*B*d
    return (xr, yr)

# --- Triangle and orthocenter reflection ----------------------------------
A, B, C = (0.0, 0.0), (2.0, 0.0), (0.5, 1.5)

# Sides as lines
BC = line_through(B,C)
CA = line_through(C,A)
AB = line_through(A,B)

# Altitudes: through A ⟂ BC, through B ⟂ CA
altA = perpendicular_through(BC, A)
altB = perpendicular_through(CA, B)

# Orthocenter H = altA ∩ altB
H = intersection(altA, altB)

# Reflect H across side BC to H'
Hprime = reflect_point_across_line(H, BC)

# --- Euclidean verification ------------------------------------------------
# is_pia_cyclic with K=0 reduces to Euclidean cyclicity check
print("Euclidean (K=0) cyclicity:", is_pia_cyclic(A,B,C, Hprime, ConstantCurvature(0.0)))

# --- Adaptive-π check (first-order) ----------------------------------------
Kb = GaussianBump(K0=0.02, x0=0.8, y0=0.4, sigma=0.45)
print("πₐ (Gaussian bump) cyclicity:", is_pia_cyclic(A,B,C, Hprime, Kb))

# Optional: numeric circumcenter consistency (Euclidean proxy)
# Compute circumcenter and radius; confirm |OH'|-R ~ 0

def circumcenter(P, Q, R):
    (ax,ay),(bx,by),(cx,cy)=P,Q,R
    d = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    if abs(d) < 1e-12:
        return ((ax+bx+cx)/3.0, (ay+by+cy)/3.0), float('nan')
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux = (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by))/d
    uy = (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax))/d
    R = math.hypot(ax-ux, ay-uy)
    return (ux,uy), R

O,R = circumcenter(A,B,C)
err = abs(math.hypot(Hprime[0]-O[0], Hprime[1]-O[1]) - R)
print("Euclidean radius error |OH'|-R:", err)
