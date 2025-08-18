#!/usr/bin/env python3
import argparse, csv
import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

def local_minima_below(arr, thr):
    arr_mid = arr[1:-1]
    is_min = (arr_mid < arr[:-2]) & (arr_mid < arr[2:]) & (arr_mid <= thr)
    idx = np.where(is_min)[0] + 1
    return idx

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--n_primes", type=int, default=100)
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--tmax", type=float, default=50.0)
    ap.add_argument("--tn", type=int, default=5000)
    ap.add_argument("--sigma_gauss", type=float, default=0.5)
    ap.add_argument("--dip_threshold", type=float, default=0.10)
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    known_t = np.array([14.13472514, 21.02203964, 25.01085758, 30.42487613, 32.93506159,
                        37.58617816, 40.91871901, 43.32707328, 48.00515088, 49.77383248])

    primes = list(primerange(2, 600))[:args.n_primes]
    logp = np.log(np.array(primes, dtype=float))
    t = np.linspace(args.tmin, args.tmax, args.tn)

    # H(t)
    H = np.zeros_like(t)
    for lp in logp: H += np.cos(2*np.pi*lp*t)
    H_norm = (H - H.min()) / (H.max() - H.min())

    # π_a line
    pi_line = np.zeros_like(t)
    for p, lp in zip(primes, logp):
        pi_line += np.exp(- (t - lp)**2 / (2*args.sigma_gauss**2)) / p
    pi_line_norm = (pi_line - pi_line.min()) / (pi_line.max() - pi_line.min())

    idx_H = local_minima_below(H_norm, args.dip_threshold)
    idx_pi = local_minima_below(1 - pi_line_norm, args.dip_threshold)

    t_H_dips = t[idx_H]; t_pi_dips = t[idx_pi]

    with open(os.path.join(args.outdir, "dip_matches.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t_H_dip","nearest_pi_a_dip","distance","nearest_zero_distance"])
        for tt in t_H_dips:
            if len(t_pi_dips) > 0:
                j = np.argmin(np.abs(t_pi_dips - tt)); d = float(np.abs(t_pi_dips[j] - tt)); tp = float(t_pi_dips[j])
            else: d = float("nan"); tp = float("nan")
            d0 = float(np.min(np.abs(known_t - tt))); w.writerow([float(tt), tp, d, d0])

    plt.figure(figsize=(10,4))
    plt.plot(t, H_norm, label="H(t) normalized")
    plt.axhline(args.dip_threshold, linestyle="--", label="dip threshold")
    for z in known_t: plt.axvline(z, linestyle=":", alpha=0.5)
    plt.scatter(t_H_dips, H_norm[idx_H], s=20, zorder=3)
    plt.title("Harmonic trace H(t) with dip candidates and known zero lines")
    plt.xlabel("t"); plt.ylabel("H_norm"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "H_trace.png"), dpi=200); plt.close()

    plt.figure(figsize=(10,4))
    plt.plot(t, pi_line_norm, label="πₐ line (normalized)")
    plt.scatter(t_pi_dips, pi_line_norm[idx_pi], s=20, zorder=3)
    plt.title("πₐ line with minima (via dips of 1 - πₐ)")
    plt.xlabel("t"); plt.ylabel("πₐ_norm"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "pi_a_line.png"), dpi=200); plt.close()

    plt.figure(figsize=(10,4))
    plt.plot(t, H_norm, label="H_norm")
    plt.plot(t, 1 - pi_line_norm, label="1 - πₐ_norm", alpha=0.7)
    plt.scatter(t_H_dips, H_norm[idx_H], s=18, label="H dips", zorder=3)
    plt.scatter(t_pi_dips, 1 - pi_line_norm[idx_pi], s=18, label="πₐ dips", zorder=3)
    for z in known_t: plt.axvline(z, linestyle=":", alpha=0.4)
    plt.title("Dip overlay: H(t) and πₐ minima with known zeros")
    plt.xlabel("t"); plt.ylabel("normalized amplitudes"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "dip_overlay.png"), dpi=220); plt.close()

if __name__ == "__main__":
    import numpy as np
    main()
