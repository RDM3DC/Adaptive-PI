# adaptive-pi-geometry

A lightweight Python library + reference spec for **Adaptive \u03c0 (\u03c0\u2090)** geometry: a unifying framework that reduces to Euclidean geometry when curvature \u2192 0 and extends gracefully to curved/inhomogeneous settings using first\u2011order Gauss\u2013Bonnet corrections.

---

## Why this repo exists

* **Goal:** Provide clean, testable primitives for \u03c0\u2090 so others can *use* and *verify* the framework on classic Euclidean (IMO-style) problems, then turn the \u201ccurvature dial\u201d to see controlled deviations.
* **Philosophy:** Everything here must (1) recover standard Euclidean results when `K(x,y) = 0`, and (2) make curvature effects explicit and composable.

---

## Repository structure

```
adaptive-pi-geometry/
├─ README.md
├─ LICENSE
├─ pyproject.toml
├─ .gitignore
├─ docs/
│  ├─ SPEC.md
│  └─ IMO-translation-guide.md
├─ src/
│  └─ pi_a/
│     ├─ __init__.py
│     ├─ core.py            # \u03c0\u2090 primitives: angle sum, cyclicity, flux integrals
│     ├─ geometry.py        # basic constructions & helpers
│     └─ models.py          # curvature fields (K) + convenience wrappers
├─ examples/
│  ├─ example1_orthocenter_reflection_pi_a.py
│  └─ example2_ceva_trig_pia.py
└─ tests/
   ├─ test_flat_limit.py
   ├─ test_cyclicity_balance.py
   └─ test_ceva_trig_pia.py
```

---

## Quickstart

```bash
# 1) clone your new repo (after you push it), then:
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .
pytest -q
```

```python
# Minimal usage
from pi_a.core import triangle_angle_sum_pia, is_pia_cyclic
from pi_a.models import ConstantCurvature

A, B, C = (0.0, 0.0), (1.0, 0.0), (0.2, 0.8)
K = ConstantCurvature(0.0)   # Flat: \u03c0\u2090 \u2192 \u03c0
S = triangle_angle_sum_pia(A,B,C, K)
print(S)  # ~ math.pi

# Check \u03c0\u2090-cyclicity for quadruple (recovers Euclidean cyclic test when K=0)
D = (0.9, 0.2)
print(is_pia_cyclic(A,B,C,D, K))
```

---

## Installation

```bash
pip install -e .
```

---

## Design guarantees

* **Flat\u2011limit exactness:** With `K(x,y) \u2261 0`, all \u03c0\u2090 primitives reduce to classical Euclidean identities.
* **First\u2011order correctness:** For small/slowly\u2011varying curvature, results match Gauss\u2013Bonnet style corrections with integrated curvature flux through polygons/sectors.
* **Composability:** Angle sums, cyclicity, Ceva/Menelaus (trig form) expose the curvature terms explicitly so you can budget/trace deviations.

---

## Files

... (See individual files for details).

---

## Roadmap

* **v0.2**: Explicit \u03c0\u2090 Ceva/Menelaus numerics (build sector/flux terms for split points).
* **v0.3**: Inversion and power\u2011of\u2011a\u2011point under \u03c0\u2090; adaptive circumcurves.
* **v0.4**: Geodesic circle model & better sector integration (adaptive radius profile).
* **v0.5**: Symbolic layer (SymPy) for first\u2011order expansions and human\u2011readable proofs.

---

## How to publish

```bash
git init
git add .
git commit -m "init: adaptive-pi-geometry v0.1 bootstrap"
git branch -M main
git remote add origin git@github.com:RDM3DC/adaptive-pi-geometry.git
git push -u origin main
```

---

## Citation (placeholder)

> Stephens, R.K. (2025). *Adaptive \u03c0 (\u03c0\u2090) Geometry: Flat\u2011Limit Exactness and First\u2011Order Curvature Corrections.* RDM3DC Labs. Version 0.1.
