
# Phase 9 — Torch-accelerated RMSE vs known zeros (t>100) + Permutation Null

This phase adds a **prediction-quality metric**: RMSE between predicted zero heights
(peak centers from the σ=0.5 closure t-density) and **known zero heights** for t>100.
We also compute a **permutation null** comparing σ=0.5 to off‑σ (e.g., 0.45, 0.55).

## Run (demo)
```bash
python scripts/phase9_sigma_torch.py --outdir outputs
```

## Heavy run
```bash
python scripts/phase9_sigma_torch.py --outdir outputs --n_primes 50000 --grid 1024   --steps 1200 --seed_grid 12 --angles 8 --tmax 140 --B 300 --use_torch 1
```

## Outputs
- `phase9_sigma05_t_hist_peaks_gt100.png` — σ=0.5 t‑density with predicted peaks + zeros (t>100)
- `phase9_summary.json` — RMSE at σ=0.5 and permutation p‑values vs off‑σ
- `phase9_tvals_sigma05.csv`, `phase9_tvals_off.csv` — raw t‑samples (for reproducibility)

Notes
- If PyTorch is installed, `--use_torch 1` accelerates base‑field building; otherwise NumPy is used.
- Increase primes/grid/steps for sharper peaks; use `--tmax 140` (or higher) to include more zeros >100.
