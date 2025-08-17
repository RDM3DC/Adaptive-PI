# Advanced πₐ Primitives

This document catalogs higher-level constructions built atop the core flux + angle primitives.

## Contents
1. Ceva & Menelaus (flux-corrected trig forms)
2. Power / Inversion (see `inversion.md`)
3. Miquel Point (flux-stabilized)
4. Nine-Point Curve (balanced circumcurve)
5. Projective Cross-Ratio (first-order invariance)
6. Cyclic Quadrilateral Opposite Angles

---
## 1. Ceva & Menelaus
Modules: `ceva.py`, `menelaus.py`

Both trig identities acquire a shared multiplicative factor exp(½ Σ ΔΦ_vertex). To first order:
  TrigCeva_πₐ = TrigCeva_Euclid * e^{Ξ},  TrigMenelaus_πₐ = TrigMenelaus_Euclid * e^{Ξ}.
With Euclidean products = 1, logical equivalence survives under curvature.

API:
- `trig_ceva_pia(A,B,C,D,E,F,K, use_flux=True)`
- `concurrent_pia(...)` boolean wrapper
- `trig_menelaus_pia(A,B,C,P,Q,R,K, use_flux=True)`
- `collinear_pia(...)`

---
## 2. Miquel Point
Module: `miquel.py`

Computes Euclidean Miquel point first, then locally refines via minimizing sum of sector flux imbalances across quadruple circles. Gradient-free finite-difference steps.

API:
- `miquel_point_euclid(A,B,C,D)`
- `miquel_point_pia(A,B,C,D,K, iters=20)` → (point, residual)

Residual near zero in flat limit; growth indicates curvature-induced drift.

---
## 3. Nine-Point Curve
Module: `ninepoint.py` + `circumcurve.py`

Construct classical 9 points (midpoints, feet, Euler midpoints). Center approximation via `circumcurve_pia_balanced` which nudges Euclidean circumcenter to reduce sector flux imbalance.

API:
- `nine_points(A,B,C)` list
- `ninepoint_curve_pia(A,B,C,K)` → (center, radius, points)

---
## 4. Projective Cross-Ratio
Module: `projective.py`

First-order claim: collinear quadruple cross-ratio is curvature invariant up to O(K^2) in normal coordinates. Current implementation leaves πₐ form identical to Euclidean; future work will add lens adjustments for large arcs.

API:
- `cross_ratio(A,B,C,D)`
- `cross_ratio_pia(A,B,C,D)` (alias for now)

---
## 5. Cyclic Quadrilateral Opposite Angles
Module: `cyclic.py`

Opposite inscribed angle sum:
  ∠ABC + ∠ADC ≈ π + ½(Φ_sector( B↔C ) - Φ_sector( D↔A )).
Boolean tester uses tolerance; future refinement will auto-estimate πₐ via triangle flux average.

API:
- `opposite_angle_sum_pia(A,B,C,D,K, center)`
- `is_pia_cyclic_quad(A,B,C,D,K, center, tol)`

---
## 6. Inversion / Power
See `inversion.md`.

---
## Shared Themes
- **Flux locality:** All corrections rely on small sectors or lenses; accuracy tied to sector sampling resolution.
- **Flat-limit sanity:** Every primitive collapses to textbook result with `K≡0`.
- **Decomposition pattern:** (Euclidean quantity) + (curvature flux) + higher-order remainder.

---
## Roadmap Enhancements
| Primitive | Next Step | Rationale |
|-----------|-----------|-----------|
| Ceva/Menelaus | Non-degenerate angle selection & variance reduction | Reduce numerical zero sines |
| Miquel | Quasi-Newton residual descent | Faster convergence |
| Nine-Point | Euler line drift quantification | Theoretical insight |
| Cross-Ratio | Add second-order correction term | Projective robustness |
| Cyclic Quad | Auto center optimization | Remove external center guess |
| Inversion | Adaptive lens triangulation | Accuracy on steep K |

---
**Status:** Experimental; API may evolve. Feedback + PRs welcome.
