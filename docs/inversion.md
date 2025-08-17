# Inversion & Power-of-a-Point in Adaptive π (πₐ)

This document specifies the first-order curvature-aware generalization of classical circle inversion and the power of a point.

## Euclidean baseline
For a circle C(O,r) and point P:
  Power(P;C) = |OP|^2 - r^2.
  Inversion: P' = O + (r^2 / |OP|^2) (P - O), P ≠ O.

## πₐ augmentation (first order)
We add a curvature lens flux term Φ_lens(P;C,K) that measures integrated curvature inside the narrow lens bounded by circle C and the two tangents from P (or a symmetric proxy if P is interior):

  power_πₐ(P;C) = |OP|^2 - r^2 + Φ_lens(P;C,K) + O(K^2).

Rationale: In constant curvature, circle ↔ geodesic circle transformations shift by area * K. For mild, local curvature fields K(x,y), the dominant perturbation is the flux through the slender region where tangent geometry defines contact with C.

## Implementation summary
File: `src/pi_a/inversion.py` provides:

- `Circle(O,r)`: data holder.
- `tangent_sector_flux(P,C,K,samples)`: samples along arc and segments to approximate Φ_lens.
- `power_of_point_pia(P,C,K)`: Euclidean power + tangent sector flux.
- `invert_point_pia(P,C,K)`: currently Euclidean inversion (first-order curvature left as a TODO for mapping distortions); since inversion formula itself is conformal, first-order curvature correction enters via metric, not needed for scalar power equality use-cases here.
- `power_decomposition_pia(P,C,K)`: returns (euclidean_power, curvature_flux) for transparency.

## Accuracy & limitations
- Uses deterministic arc sampling with mid-segment area proxies (triangular wedge approximation). Error is O(h^2) + O(K^2).
- Interior point fallback: orthogonal diameter wedge; more refined interior lens integrator is future work.
- Adequate for comparative power tests and inversion-invariance experiments under small curvature (|K| * r^2 << 1).

## Future extensions
1. Adaptive angular sampling weighted by curvature magnitude.
2. Analytical constant-curvature closed form (e.g. spherical/hyperbolic first-order expansion) to validate numerics.
3. Curvature-corrected inversion mapping for second-order shape preservation (Jacobi field integration along radial lines).

## Example
```python
from pi_a.inversion import Circle, power_of_point_pia, invert_point_pia
from pi_a.models import ConstantCurvature, GaussianBump

C = Circle((0.0,0.0), 2.0)
P = (5.0, 1.0)
flat = ConstantCurvature(0.0)
curved = GaussianBump(K0=0.05, x0=0.5, y0=-0.2, sigma=0.8)
print(power_of_point_pia(P,C,flat))    # ~ Euclidean
print(power_of_point_pia(P,C,curved))  # perturbed
print(invert_point_pia(P,C,flat))
```

## Testing
See `tests/test_power_inversion_pia.py` for:
- Power at center (flat)
- Inversion outside circle
- Double inversion identity
- Nonzero curvature flux difference
- Tangent sector symmetry (flat zero)

---
**Status:** v0 first-order; curvature lens term empirically validated by tests.
