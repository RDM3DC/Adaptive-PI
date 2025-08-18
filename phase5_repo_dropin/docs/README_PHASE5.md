
# Phase 5 — πₐ RH Prototypes

This folder contains two scripts for the Riemann–πₐ experiments.

## 1) `scripts/phase5_geodesic_tracers.py`

Integrates geodesic-like tracers over a 2D πₐ field (prime-weighted Gaussians) and detects near-closures. Each closure's average y-position maps to an imaginary-part estimate `t` which is compared to the first 10 known nontrivial ζ zeros.

**Run**:
```bash
python scripts/phase5_geodesic_tracers.py --outdir outputs
# knobs:
#   --grid 280 --primes 200 --steps 1100 --closure_tol 0.03 --smooth_px 2.5
```

**Outputs**:
- `pia_field_trajectories.png` — φ field + sample trajectories
- `closure_t_hist.png` — histogram of t estimates with known zero lines
- `closure_matches.csv` — closure records with nearest-zero distance

## 2) `scripts/phase5_harmonic_overlay.py`

Builds a 1D harmonic trace
\[ H(t)=\sum\limits_{p\le N} \cos(2\pi \log p\, t) \]
and a πₐ line via Gaussian bumps at `log p`, then overlays their dip locations with the known zero lines.

**Run**:
```bash
python scripts/phase5_harmonic_overlay.py --outdir outputs
# knobs:
#   --n_primes 200 --tn 8000 --sigma_gauss 0.35 --dip_threshold 0.08
```

**Outputs**:
- `H_trace.png`, `pi_a_line.png`, `dip_overlay.png`, `dip_matches.csv`

---

## Tips
- More primes + larger grids sharpen structure (use with care for runtime).
- For formal runs, log your parameters and seed grids in the CSVs.
- Compare Phase 5 results against Phase II.5 σ-sweep closure dips for a 2-pronged validation.
