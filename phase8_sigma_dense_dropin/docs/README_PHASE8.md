
# Phase 8 — Dense σ‑sweep + Zero‑prediction metrics

This step scales the conformal‑geodesic program and adds quantitative **zero‑prediction** tests.

## Run (demo)
```bash
python scripts/phase8_sigma_dense.py --outdir outputs
```

## Heavy run (as per X suggestions)
```bash
python scripts/phase8_sigma_dense.py --outdir outputs \
  --n_primes 10000 --grid 1024 \
  --sigma_min 0.25 --sigma_max 0.75 --sigma_step 0.025 \
  --B 300 --steps 1200 --smooth_px 2.5 --seed_grid 12 --angles 8 --tmax 100
```

### Outputs
- `phase8_mean_closure_ci.png` — mean closure error vs σ with 95% bootstrap bands
- `phase8_count_ci.png` — closure count vs σ with 95% bootstrap bands
- `phase8_closure_density_heatmap.png` — closure density over (σ, t)
- `phase8_sigma_results.json`, `phase8_closures.csv` — raw sweep data
- `phase8_zero_metrics.json` — zero‑prediction metrics (accuracy@τ, distance stats, peak hits)
- `phase8_t_hist_sigma05_peaks.png` — σ=0.5 histogram with known zeros and detected density peaks

Date: 2025-08-18
