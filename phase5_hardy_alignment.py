#!/usr/bin/env python3
"""
Phase 5 — Hardy-style alignment:
- δ-smoothed π_a(x) = π + Σ (1/p) δ_σ(x - log p), implemented with Gaussian smoothing
- φ(x) = -log π_a(x); curvature proxy K ≈ -φ''
- H_{1/2}(t) = Re Σ p^{-1/2} e^{i t log p} vs Hardy Z(t) on the critical line
- Detect dips of H (bottom q-th percentile), detect Z zeros (sign changes)
- Cross-correlate H and Z; export plots + CSV

Usage:
  python scripts/phase5_hardy_alignment.py --outdir outputs
Knobs:
  --n_primes 400 --sigma 0.1 --xmax 7.0 --nx 3000 --tmax 80 --nt 3000 --dip_q 0.10
"""
import argparse, csv, os
import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange
import mpmath as mp

def build_pi_a(logp, x, sigma):
    pi_a = np.full_like(x, np.pi, dtype=float)
    inv2s2 = 1.0 / (2 * sigma * sigma)
    for p, lp in logp:
        pi_a += (1.0 / p) * np.exp(-(x - lp) ** 2 * inv2s2)
    return pi_a

def hardy_theta(tval):
    z = 0.25 + 0.5j * tval
    return float(mp.im(mp.loggamma(z)) - 0.5 * tval * mp.log(mp.pi))

def hardy_Z(tval):
    s = 0.5 + 1j * tval
    theta = hardy_theta(tval)
    return float(mp.re(mp.e**(1j*theta) * mp.zeta(s)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--n_primes", type=int, default=200)
    ap.add_argument("--sigma", type=float, default=0.15, help="Gaussian width for δ smoothing on x=log p")
    ap.add_argument("--xmin", type=float, default=0.0)
    ap.add_argument("--xmax", type=float, default=6.0)
    ap.add_argument("--nx", type=int, default=2000)
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--tmax", type=float, default=50.0)
    ap.add_argument("--nt", type=int, default=1500)
    ap.add_argument("--dip_q", type=float, default=0.10, help="percentile for H dips (0.10 = bottom 10%)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    mp.mp.dps = 50

    # Primes and logs
    # Upper bound heuristic for primerange so we have enough primes to slice
    upper = max(1000, int(np.exp(args.xmax)) + 50)
    primes = list(primerange(2, upper))[:args.n_primes]
    logp = [(int(p), float(np.log(p))) for p in primes]

    # π_a, φ, K proxy over x ~ log p axis
    x = np.linspace(args.xmin, args.xmax, args.nx)
    pi_a = build_pi_a(logp, x, args.sigma)
    phi = -np.log(pi_a)
    dx = (args.xmax - args.xmin) / (args.nx - 1)
    phi_prime = np.gradient(phi, dx)
    phi_second = np.gradient(phi_prime, dx)
    K_proxy = -phi_second

    # H_{1/2}(t) & Hardy Z(t)
    t = np.linspace(args.tmin, args.tmax, args.nt)
    H = np.zeros_like(t)
    for p, lp in logp:
        H += (p ** -0.5) * np.cos(t * lp)
    H_z = (H - H.mean()) / (H.std() + 1e-12)

    Z = np.array([hardy_Z(tv) for tv in t], dtype=float)
    Z_z = (Z - Z.mean()) / (Z.std() + 1e-12)

    # Dips and zeros
    thr = np.quantile(H_z, args.dip_q)
    is_min = np.r_[False, (H_z[1:-1] < H_z[:-2]) & (H_z[1:-1] < H_z[2:]), False]
    dip_idx = np.where(is_min & (H_z <= thr))[0]
    t_H_dips = t[dip_idx]

    signs = np.sign(Z_z)
    zcross = np.where(signs[:-1] * signs[1:] < 0)[0]
    t_Z_zeros = 0.5 * (t[zcross] + t[zcross + 1]) if len(zcross) else np.array([])

    # Cross-correlation
    Hn = H_z - H_z.mean()
    Zn = Z_z - Z_z.mean()
    corr = np.correlate(Hn, Zn, mode='full')
    lags = np.arange(-len(Hn) + 1, len(Hn))
    imax = int(np.argmax(corr))
    best_lag = int(lags[imax])
    best_corr = float(corr[imax] / (np.linalg.norm(Hn) * np.linalg.norm(Zn) + 1e-12))

    # CSV of dip→nearest zero distances
    csv_path = os.path.join(args.outdir, "hardy_alignment_dips_vs_zeros.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_H_dip", "nearest_Z_zero", "abs_diff"])
        for th in t_H_dips:
            if len(t_Z_zeros):
                j = int(np.argmin(np.abs(t_Z_zeros - th)))
                tz = float(t_Z_zeros[j])
                w.writerow([float(th), tz, float(abs(tz - th))])

    # Plots
    # 1) π_a, φ, K proxy
    plt.figure(figsize=(9,4))
    plt.plot(x, pi_a, label="pi_a(x)")
    plt.plot(x, phi, label="phi(x) = -log pi_a")
    plt.plot(x, K_proxy, label="K proxy ≈ -phi''")
    plt.title("δ-smoothed π_a, potential φ, and curvature proxy (K ≈ -φ'')")
    plt.xlabel("x  (≈ log p)")
    plt.ylabel("value")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "pia_phi_K_profile.png"), dpi=220)
    plt.close()

    # 2) Overlay H_{1/2}(t) and Hardy Z(t)
    plt.figure(figsize=(10,4))
    plt.plot(t, H_z, label="H_{1/2}(t) (normalized)")
    plt.plot(t, Z_z, label="Hardy Z(t) (normalized)")
    for tz in t_Z_zeros:
        plt.axvline(tz, linestyle=":", alpha=0.35)
    plt.title("H_{1/2}(t) vs Hardy Z(t) (normalized) with zero crossings")
    plt.xlabel("t")
    plt.ylabel("amplitude (z-scored)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "H_vs_Z_overlay.png"), dpi=220)
    plt.close()

    # 3) Cross-correlation curve
    plt.figure(figsize=(9,4))
    denom = (np.linalg.norm(Hn) * np.linalg.norm(Zn) + 1e-12)
    plt.plot(lags, corr / denom)
    plt.title(f"Cross-correlation H vs Z  (max at lag={best_lag}; r≈{best_corr:.3f})")
    plt.xlabel("lag (samples)")
    plt.ylabel("normalized correlation")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "HZ_crosscorr.png"), dpi=220)
    plt.close()

    # 4) Mark dips and zeros on H trace
    plt.figure(figsize=(10,4))
    plt.plot(t, H_z, label="H_{1/2}(t) (normalized)")
    if len(dip_idx):
        plt.scatter(t_H_dips, H_z[dip_idx], s=12, label="H dips")
    for tz in t_Z_zeros:
        plt.axvline(tz, linestyle=":", alpha=0.35)
    plt.title(f"H dips (bottom {int(args.dip_q*100)}%) vs Z(t) zeros")
    plt.xlabel("t")
    plt.ylabel("H z-score")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "H_dips_vs_Z_zeros.png"), dpi=220)
    plt.close()

    # Print a tiny summary
    print(f"Max corr r≈{best_corr:.3f} at lag {best_lag}; wrote CSV to {csv_path}")

if __name__ == "__main__":
    main()
