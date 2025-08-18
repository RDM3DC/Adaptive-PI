#!/usr/bin/env python3
"""
Phase 8 — Dense σ-sweep with bootstrap CIs and zero‑prediction metrics.

Features
- Conformal geodesics on π_a field (prime-weighted Gaussians)
- Dense σ-sweep (range/step configurable)
- Bootstrap bands for mean closure error & closure counts
- (σ, t) closure heatmap
- Zero-prediction metrics at σ=0.5:
  * accuracy@τ (fraction of closures within τ of any known zero)
  * distance quantiles (mean/median/P25/P75)
  * peak detection on a smoothed t-density and hit‑rate vs known zeros

Heavy runs: n_primes=10000, grid_n=1024, steps≈1200, B≈300 (use your workstation).

Usage (light demo):
  python scripts/phase8_sigma_dense.py --outdir outputs

Suggested heavy run:
  python scripts/phase8_sigma_dense.py --outdir outputs --n_primes 10000 --grid 1024 \
      --sigma_min 0.25 --sigma_max 0.75 --sigma_step 0.025 --B 300 --steps 1200 --smooth_px 2.5
"""
import argparse, os, csv, json, math
import numpy as np
import matplotlib.pyplot as plt
from math import log

# ---------- Utilities ----------
def sieve_primes(n: int):
    sieve = np.ones(n+1, dtype=bool); sieve[:2] = False
    for p in range(2, int(n**0.5)+1):
        if sieve[p]: sieve[p*p:n+1:p] = False
    return [int(i) for i in range(n+1) if sieve[i]]

def gaussian_kernel1d(sigma_px: float):
    radius = max(1, int(3*sigma_px))
    xs = np.arange(-radius, radius+1)
    k = np.exp(-0.5*(xs/sigma_px)**2); k /= k.sum()
    return k

def smooth_gaussian(field: np.ndarray, sigma_px: float = 2.0):
    ker = gaussian_kernel1d(sigma_px)
    tmp = np.apply_along_axis(lambda m: np.convolve(m, ker, mode="same"), axis=1, arr=field)
    out = np.apply_along_axis(lambda m: np.convolve(m, ker, mode="same"), axis=0, arr=tmp)
    return out

def bilinear_sample(field: np.ndarray, x: float, y: float) -> float:
    h, w = field.shape; eps = 1e-9
    x = max(0.0, min(1.0 - eps, x)); y = max(0.0, min(1.0 - eps, y))
    gx = x * (w - 1); gy = y * (h - 1)
    x0 = int(np.floor(gx)); y0 = int(np.floor(gy))
    x1 = min(x0 + 1, w - 1); y1 = min(y0 + 1, h - 1)
    tx = gx - x0; ty = gy - y0
    v00 = field[y0, x0]; v10 = field[y0, x1]; v01 = field[y1, x0]; v11 = field[y1, x1]
    return (1-ty)*((1-tx)*v00 + tx*v10) + ty*((1-tx)*v01 + tx*v11)

# ---------- π_a base field ----------
def build_base_prime_field(grid_n: int = 300, n_primes: int = 600, bump_sigma: float = 0.028):
    # upper bound heuristic for sieve
    sieve_n = max(1000, int(1.3 * n_primes * np.log(max(10, n_primes))))
    primes = sieve_primes(sieve_n)[:n_primes]
    pmax = primes[-1]; logpmax = log(pmax)
    ys, xs = np.mgrid[0:grid_n, 0:grid_n]
    X = xs / (grid_n - 1); Y = ys / (grid_n - 1)
    inv2s2 = 1.0 / (2.0 * bump_sigma * bump_sigma)
    base = np.zeros((grid_n, grid_n), dtype=np.float64)
    for p in primes:
        xp = p / pmax; yp = log(p) / logpmax
        base += np.exp(-((X - xp)**2 + (Y - yp)**2) * inv2s2) / p
    return base

# ---------- Conformal geodesic integrator ----------
def geodesic_rk4(x0: float, y0: float, theta: float, phi_gx: np.ndarray, phi_gy: np.ndarray,
                 dt: float = 0.004, steps: int = 800):
    v = np.array([np.cos(theta), np.sin(theta)], dtype=np.float64)
    pos = np.array([x0, y0], dtype=np.float64)
    path = [pos.copy()]

    def grad_phi_at(px, py):
        return np.array([bilinear_sample(phi_gx, px, py), bilinear_sample(phi_gy, px, py)], dtype=np.float64)

    for _ in range(steps):
        def acc(p, v):
            g = grad_phi_at(p[0], p[1])
            vv = float(np.dot(v, v))
            return -2.0 * np.dot(g, v) * v + g * vv

        p1 = pos; v1 = v; a1 = acc(p1, v1)
        p2 = pos + 0.5*dt*v1; v2 = v + 0.5*dt*a1; a2 = acc(p2, v2)
        p3 = pos + 0.5*dt*v2; v3 = v + 0.5*dt*a2; a3 = acc(p3, v3)
        p4 = pos + dt*v3;     v4 = v + dt*a3;     a4 = acc(p4, v4)

        pos = pos + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
        v   = v   + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

        if pos[0] < 0.0 or pos[0] > 1.0:
            v[0] *= -1.0; pos[0] = float(np.clip(pos[0], 0.0, 1.0))
        if pos[1] < 0.0 or pos[1] > 1.0:
            v[1] *= -1.0; pos[1] = float(np.clip(pos[1], 0.0, 1.0))
        path.append(pos.copy())
    return np.array(path)

def closure_error(path: np.ndarray, warmup: int = 160):
    if len(path) <= warmup: return 1e9, warmup
    start = path[0]; seg = path[warmup:]
    d = np.linalg.norm(seg - start, axis=1); i = int(np.argmin(d))
    return float(d[i]), warmup + i

def build_phi_from_base(base_field: np.ndarray, sigma: float, eps: float, alpha: float, smooth_px: float):
    bias = 1.0 / (abs(float(sigma) - 0.5) + eps)
    pi_a = np.pi + bias * base_field
    pi_a = smooth_gaussian(pi_a, sigma_px=smooth_px)
    centered = pi_a - np.pi
    centered /= (np.max(np.abs(centered)) + 1e-12)
    return alpha * centered

# ---------- Known zeros (first 50 imag parts) ----------
KNOWN_T = np.array([
  14.13472514, 21.02203964, 25.01085758, 30.42487613, 32.93506159,
  37.58617816, 40.91871901, 43.32707328, 48.00515088, 49.77383248,
  52.97032148, 56.44624770, 59.34704400, 60.83177852, 65.11254405,
  67.07981053, 69.54640171, 72.06715767, 75.70469070, 77.14484007,
  79.33737502, 82.91038085, 84.73549298, 87.42527461, 88.80911121,
  92.49189927, 94.65134404, 95.87063423, 98.83119422, 101.3178510,
  103.7255380, 105.4466230, 107.1686112, 111.0295355, 111.8746592,
  114.3202209, 116.2266803, 118.7907829, 121.3701250, 122.9468293,
  124.2568186, 127.5166839, 129.5787042, 131.0876885, 133.4977372,
  134.7565097, 138.1160421, 139.7362089, 141.1237074, 143.1118458
])

# ---------- Zero-prediction metrics ----------
def accuracy_at_tau(t_vals, zeros, tau):
    if len(t_vals) == 0: return 0.0
    d = np.min(np.abs(t_vals[:,None] - zeros[None,:]), axis=1)
    return float(np.mean(d <= tau))

def dist_stats(t_vals, zeros):
    if len(t_vals) == 0:
        return {"mean": None, "median": None, "p25": None, "p75": None}
    d = np.min(np.abs(t_vals[:,None] - zeros[None,:]), axis=1)
    return {
        "mean": float(np.mean(d)),
        "median": float(np.median(d)),
        "p25": float(np.percentile(d, 25)),
        "p75": float(np.percentile(d, 75))
    }

def peak_hits(t_vals, zeros, T_MAX, bins=200, smooth_px=2.0, k_top=20, tau=0.25):
    if len(t_vals) == 0:
        return {"hits": 0, "total": 0, "rate": 0.0, "peaks": []}
    hist, edges = np.histogram(t_vals, bins=bins, range=(0.0, T_MAX))
    ker = gaussian_kernel1d(smooth_px*1.5)
    sm = np.convolve(hist, ker, mode="same")
    # local maxima
    idx = np.where((sm[1:-1] > sm[:-2]) & (sm[1:-1] > sm[2:]))[0] + 1
    # pick top-k
    order = np.argsort(sm[idx])[::-1]
    idx = idx[order][:k_top]
    centers = 0.5*(edges[idx] + edges[idx+1])
    # hit if within tau of any zero
    hits = 0
    for c in centers:
        if np.min(np.abs(zeros - c)) <= tau:
            hits += 1
    return {"hits": int(hits), "total": int(len(centers)), "rate": float(hits/max(1,len(centers))), "peaks": centers.tolist()}

# ---------- Bootstrap σ-sweep ----------
def sigma_bootstrap_sweep(base_field: np.ndarray, sigmas, seeds, angles,
                          B: int, q_closure: float, steps: int, dt: float,
                          eps: float, alpha: float, smooth_px: float, T_MAX: float,
                          outdir: str):
    os.makedirs(outdir, exist_ok=True)
    results = {}
    all_closures = []

    seed_list = [(x0, y0, th) for (x0, y0) in seeds for th in angles]
    N = len(seed_list)

    for s in sigmas:
        phi = build_phi_from_base(base_field, s, eps, alpha, smooth_px)
        gy, gx = np.gradient(phi)

        errs = np.zeros(N, dtype=float)
        t_ests = np.zeros(N, dtype=float)
        for i, (x0, y0, th) in enumerate(seed_list):
            path = geodesic_rk4(x0, y0, th, gx, gy, dt=dt, steps=steps)
            err, idx = closure_error(path, warmup=160)
            errs[i] = err
            win = path[max(0, idx-25):idx+1, 1]
            t_ests[i] = float(np.mean(win) * T_MAX)

        thr = float(np.quantile(errs, q_closure))
        mask = errs <= thr
        base_count = int(mask.sum())

        # Bootstrap
        rng = np.random.default_rng(12345)
        mean_err_bs = []; count_bs = []
        for _ in range(B):
            ii = rng.integers(0, N, size=N)
            e = errs[ii]
            thr_b = float(np.quantile(e, q_closure))
            mean_err_bs.append(float(e.mean()))
            count_bs.append(int(np.sum(e <= thr_b)))
        lo_m, hi_m = np.percentile(mean_err_bs, [2.5, 97.5])
        lo_c, hi_c = np.percentile(count_bs, [2.5, 97.5])

        # Store
        results[float(s)] = {
            "mean_err": float(errs.mean()),
            "mean_err_ci": [float(lo_m), float(hi_m)],
            "count": base_count,
            "count_ci": [float(lo_c), float(hi_c)],
            "q_threshold": thr
        }

        # Persist closures for σ, and compute zero metrics at σ=0.5 later
        for (x0, y0, th), e, t in zip(seed_list, errs, t_ests):
            if e <= thr:
                all_closures.append((float(s), float(x0), float(y0), float(th), float(e), float(t)))

    # Save sweep
    with open(os.path.join(outdir, "phase8_sigma_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    with open(os.path.join(outdir, "phase8_closures.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sigma","x0","y0","theta","closure_err","t_est"])
        for row in all_closures: w.writerow(list(row))

    return results, all_closures

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--n_primes", type=int, default=500)
    ap.add_argument("--grid", type=int, default=320)
    ap.add_argument("--bump_sigma", type=float, default=0.028)
    ap.add_argument("--sigma_min", type=float, default=0.30)
    ap.add_argument("--sigma_max", type=float, default=0.70)
    ap.add_argument("--sigma_step", type=float, default=0.05)
    ap.add_argument("--B", type=int, default=120)
    ap.add_argument("--q_closure", type=float, default=0.12)
    ap.add_argument("--steps", type=int, default=700)
    ap.add_argument("--dt", type=float, default=0.004)
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--alpha", type=float, default=0.35)
    ap.add_argument("--smooth_px", type=float, default=2.2)
    ap.add_argument("--tmax", type=float, default=80.0)
    ap.add_argument("--seed_grid", type=int, default=8, help="seed grid per side (seed_grid^2 seeds)")
    ap.add_argument("--angles", type=int, default=6, help="number of launch angles")
    ap.add_argument("--tau", type=float, default=0.25, help="tolerance for zero hit checks")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Build base field
    base = build_base_prime_field(grid_n=args.grid, n_primes=args.n_primes, bump_sigma=args.bump_sigma)

    # Seeds & angles
    xs = np.linspace(0.12, 0.88, args.seed_grid); ys = np.linspace(0.12, 0.88, args.seed_grid)
    seeds = [(float(x), float(y)) for x in xs for y in ys]
    angles = [k * 2*np.pi/args.angles for k in range(args.angles)]

    # σ grid
    sigmas = np.arange(args.sigma_min, args.sigma_max + 1e-12, args.sigma_step)

    # Run sweep with bootstrap
    results, closures = sigma_bootstrap_sweep(
        base_field=base, sigmas=sigmas, seeds=seeds, angles=angles,
        B=args.B, q_closure=args.q_closure, steps=args.steps, dt=args.dt,
        eps=args.eps, alpha=args.alpha, smooth_px=args.smooth_px, T_MAX=args.tmax,
        outdir=args.outdir
    )

    # --- Plots: mean error & count with CIs ---
    xsig = sorted(results.keys())
    mean_err = [results[s]["mean_err"] for s in xsig]
    lo_m = [results[s]["mean_err_ci"][0] for s in xsig]
    hi_m = [results[s]["mean_err_ci"][1] for s in xsig]
    count_m = [results[s]["count"] for s in xsig]
    lo_c = [results[s]["count_ci"][0] for s in xsig]
    hi_c = [results[s]["count_ci"][1] for s in xsig]

    plt.figure(figsize=(8,4))
    plt.plot(xsig, mean_err, marker="o")
    plt.fill_between(xsig, lo_m, hi_m, alpha=0.2)
    plt.axvline(0.5, linestyle="--")
    plt.xlabel("σ = Re(s)"); plt.ylabel("Mean closure error (95% CI)")
    plt.title("Phase 8: Mean closure error vs σ (bootstrap)"); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "phase8_mean_closure_ci.png"), dpi=220); plt.close()

    plt.figure(figsize=(8,4))
    plt.plot(xsig, count_m, marker="o")
    plt.fill_between(xsig, lo_c, hi_c, alpha=0.2)
    plt.axvline(0.5, linestyle="--")
    plt.xlabel("σ = Re(s)"); plt.ylabel("Closure count (95% CI)")
    plt.title("Phase 8: Closure count vs σ (bootstrap)"); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "phase8_count_ci.png"), dpi=220); plt.close()

    # --- Heatmap over (σ, t) ---
    t_bins = np.linspace(0.0, args.tmax, 100)
    s_vals = np.array(xsig)
    H = np.zeros((len(t_bins)-1, len(s_vals)))
    for s, x0, y0, th, e, t in closures:
        i = int(np.argmin(np.abs(s_vals - s)))
        k = np.searchsorted(t_bins, t) - 1
        if 0 <= k < len(t_bins)-1: H[k,i] += 1.0
    plt.figure(figsize=(8,6))
    extent = [s_vals.min(), s_vals.max(), t_bins[0], t_bins[-1]]
    plt.imshow(H, origin='lower', extent=extent, aspect='auto')
    plt.colorbar(label="closure density")
    plt.axvline(0.5, linestyle="--")
    plt.xlabel("σ = Re(s)"); plt.ylabel("t"); plt.title("Phase 8: Closure density over (σ, t)")
    plt.tight_layout(); plt.savefig(os.path.join(args.outdir, "phase8_closure_density_heatmap.png"), dpi=220); plt.close()

    # --- Zero-prediction metrics at σ=0.5 ---
    t05 = np.array([t for (s, x0, y0, th, e, t) in closures if abs(s - 0.5) < 1e-9])
    zeros = KNOWN_T[KNOWN_T <= args.tmax]
    stats = dist_stats(t05, zeros)
    acc_tau = {tau: accuracy_at_tau(t05, zeros, tau) for tau in [0.1, 0.2, 0.25, 0.3, 0.5]}
    peaks = peak_hits(t05, zeros, args.tmax, bins=240, smooth_px=2.0, k_top=25, tau=args.tau)

    with open(os.path.join(args.outdir, "phase8_zero_metrics.json"), "w") as f:
        json.dump({"dist_stats": stats, "accuracy_at_tau": acc_tau, "peak_hits": peaks}, f, indent=2)

    # Plot t05 histogram + zeros + detected peaks
    plt.figure(figsize=(9,4))
    if len(t05) > 0:
        plt.hist(t05, bins=60, alpha=0.85)
    for z in zeros: plt.axvline(z, linestyle=":", alpha=0.5)
    if len(peaks["peaks"]) > 0:
        for c in peaks["peaks"]:
            plt.axvline(c, linestyle="--", alpha=0.6)
    plt.xlabel("t at σ=0.5"); plt.ylabel("count"); plt.title("Phase 8: t‑density (σ=0.5) with zeros and detected peaks")
    plt.tight_layout(); plt.savefig(os.path.join(args.outdir, "phase8_t_hist_sigma05_peaks.png"), dpi=220); plt.close()

if __name__ == "__main__":
    main()
