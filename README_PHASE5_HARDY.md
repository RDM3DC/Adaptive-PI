
# Phase 5 — Hardy Alignment Prototype

This script tests the rigorous-style sketch:
- Define a δ-smoothed \(\pi_a(x)=\pi + \sum_p \frac{1}{p}\,\delta_\sigma(x-\log p)\).
- Set \(\phi(x)=-\log\pi_a(x)\); curvature proxy \(K\approx-\phi''(x)\).
- Compare \(H_{1/2}(t)=\Re\sum p^{-1/2} e^{it\log p}\) with **Hardy** \(Z(t)\).

## Run
```bash
python scripts/phase5_hardy_alignment.py --outdir outputs
# knobs:
#   --n_primes 400 --sigma 0.1 --xmax 7.0 --nx 3000 --tmax 80 --nt 3000 --dip_q 0.10
```

## Outputs
- `pia_phi_K_profile.png` — δ-smoothed \(\pi_a\), \(\phi\), and curvature proxy
- `H_vs_Z_overlay.png` — overlay of \(H_{1/2}(t)\) and Hardy \(Z(t)\) with zero crossings
- `HZ_crosscorr.png` — cross-correlation curve
- `H_dips_vs_Z_zeros.png` — H dips (bottom q%) vs Z zeros
- `hardy_alignment_dips_vs_zeros.csv` — dip-to-zero distances

## Notes
- Increase primes and extend `tmax` for stronger alignment (cost ↑).
- Tune `sigma` (δ-width) to balance smoothing vs fidelity.
- Next: lift to 2D \(\phi(\sigma,t)\) and integrate conformal geodesics (see Phase II.5).
