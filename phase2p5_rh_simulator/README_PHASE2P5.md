# Phase II.5 — Adaptive π (πₐ) RH Geodesic Simulator

This script reproduces the figures we generated in chat:

- **closure_vs_sigma.png** — average geodesic closure error across σ = Re(s)
- **field_sigma05.png** — φ field snapshot at σ = 0.5
- **geodesics_sigma05.png** — geodesic trajectories overlay on φ
- **closure_heatmap_sigma05.png** — closure error heatmap (σ = 0.5)

## How it works

- Builds a weighted πₐ field:

  πₐ = π + Σ_p exp(-||(x - p_norm, y - log p_norm)||² / (2σ²)) / (p * |σ - 0.5| + ε)

- Smooths it with a Gaussian (no SciPy dependency) and centers to φ.
- Treats e^{2φ} δ_ij as the metric; integrates geodesics using the conformal Christoffel form.

- Measures closure error as min distance back to start after a warm-up period.


## Run

```bash
python scripts/run_phase2p5.py --outdir outputs
```

Outputs will be saved to `outputs/`.

## Notes

- Increase `--primes` and internal grid to strengthen the σ=0.5 dip.
- For reproducibility, the prime placement is deterministic.
- This is *experimental* and intended as a research prototype for πₐ-geometry × RH.
