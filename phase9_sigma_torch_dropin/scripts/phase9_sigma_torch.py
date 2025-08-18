#!/usr/bin/env python3
"""
Phase 9 — Torch-accel (optional) RMSE vs zeros (t>100) + off-line permutation null.

Key additions
- Focus on σ=0.5 predictions for zeros with imag part t>100.
- Predict zero heights from t-density peaks of near-closures at σ=0.5.
- Compute RMSE to nearest known zeros (t>100).
- Build permutation null: compare σ=0.5 RMSE to off-σ (e.g., 0.45, 0.55) by shuffling labels.
- Optional PyTorch acceleration for field building and batched geodesics.

Usage (demo):
  python scripts/phase9_sigma_torch.py --outdir outputs

Suggested heavy run:
  python scripts/phase9_sigma_torch.py --outdir outputs --n_primes 50000 --grid 1024 \
    --steps 1200 --seed_grid 12 --angles 8 --tmax 140 --B 300 --use_torch 1
"""
import argparse, os, json, csv, math, numpy as np
import matplotlib.pyplot as plt
from math import log

try:
    import torch
    TORCH_OK = True
except Exception:
    TORCH_OK = False

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
def build_base_prime_field_numpy(grid_n: int, n_primes: int, bump_sigma: float):
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

def build_base_prime_field_torch(grid_n: int, n_primes: int, bump_sigma: float, device="cpu"):
    sieve_n = max(1000, int(1.3 * n_primes * np.log(max(10, n_primes))))
    primes = sieve_primes(sieve_n)[:n_primes]
    pmax = primes[-1]; logpmax = math.log(pmax)
    xs = torch.linspace(0, 1, grid_n, device=device)
    ys = torch.linspace(0, 1, grid_n, device=device)
    X, Y = torch.meshgrid(xs, ys, indexing="xy")
    inv2s2 = 1.0 / (2.0 * bump_sigma * bump_sigma)
    base = torch.zeros((grid_n, grid_n), dtype=torch.float64, device=device)
    # Chunk primes to control memory
    B = 2000
    for i in range(0, len(primes), B):
        chunk = primes[i:i+B]
        p = torch.tensor(chunk, dtype=torch.float64, device=device)
        xp = p / pmax
        yp = torch.log(p) / logpmax
        # Broadcast
        Xb = X[...,None]; Yb = Y[...,None]
        w = torch.exp(-((Xb - xp)**2 + (Yb - yp)**2) * inv2s2) / p
        base += w.sum(dim=-1)
    return base.cpu().numpy()

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

# ---------- Known zeros (first 100 imag parts) ----------
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

def nearest_zero_rmse(pred, zeros):
    if len(pred)==0 or len(zeros)==0: return None
    d = np.min(np.abs(pred[:,None] - zeros[None,:]), axis=1)
    return float(np.sqrt(np.mean(d**2)))

def predict_peaks(t_vals, T_MAX, bins=240, smooth_px=2.0, k_top=25, t_min=100.0):
    if len(t_vals)==0: return np.array([])
    hist, edges = np.histogram(t_vals, bins=bins, range=(0.0, T_MAX))
    ker = gaussian_kernel1d(smooth_px*1.5)
    sm = np.convolve(hist, ker, mode="same")
    idx = np.where((sm[1:-1] > sm[:-2]) & (sm[1:-1] > sm[2:]))[0] + 1
    order = np.argsort(sm[idx])[::-1]
    idx = idx[order][:k_top]
    centers = 0.5*(edges[idx] + edges[idx+1])
    return centers[centers >= t_min]

def permutation_null_rmse(t_on, t_off, zeros, R=500, T_MAX=140.0, t_min=100.0):
    # RMSE difference: on - off (more negative => better on)
    pred_on = predict_peaks(t_on, T_MAX, t_min=t_min); pred_off = predict_peaks(t_off, T_MAX, t_min=t_min)
    rmse_on = nearest_zero_rmse(pred_on, zeros); rmse_off = nearest_zero_rmse(pred_off, zeros)
    if (rmse_on is None) or (rmse_off is None):
        return {"rmse_on": rmse_on, "rmse_off": rmse_off, "p_perm": None}
    obs = rmse_on - rmse_off
    # Permute labels by shuffling pooled t-values and splitting by original sizes
    pool = np.concatenate([t_on, t_off]); n_on = len(t_on); n_off = len(t_off)
    rng = np.random.default_rng(1234)
    diffs = []
    for _ in range(R):
        rng.shuffle(pool)
        a = pool[:n_on]; b = pool[n_on:n_on+n_off]
        pa = predict_peaks(a, T_MAX, t_min=t_min); pb = predict_peaks(b, T_MAX, t_min=t_min)
        ra = nearest_zero_rmse(pa, zeros); rb = nearest_zero_rmse(pb, zeros)
        if (ra is None) or (rb is None):
            continue
        diffs.append(ra - rb)
    if len(diffs)==0:
        p = None
    else:
        diffs = np.array(diffs)
        p = float((np.sum(np.abs(diffs) >= abs(obs)) + 1) / (len(diffs) + 1))
    return {"rmse_on": rmse_on, "rmse_off": rmse_off, "p_perm": p}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--n_primes", type=int, default=800)
    ap.add_argument("--grid", type=int, default=300)
    ap.add_argument("--bump_sigma", type=float, default=0.028)
    ap.add_argument("--steps", type=int, default=700)
    ap.add_argument("--dt", type=float, default=0.004)
    ap.add_argument("--smooth_px", type=float, default=2.2)
    ap.add_argument("--tmax", type=float, default=140.0)
    ap.add_argument("--seed_grid", type=int, default=8)
    ap.add_argument("--angles", type=int, default=6)
    ap.add_argument("--q_closure", type=float, default=0.12)
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--alpha", type=float, default=0.35)
    ap.add_argument("--use_torch", type=int, default=0)
    ap.add_argument("--off_sigmas", type=float, nargs="+", default=[0.45, 0.55, 0.40, 0.60])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Base field
    if args.use_torch and TORCH_OK:
        base = build_base_prime_field_torch(grid_n=args.grid, n_primes=args.n_primes, bump_sigma=args.bump_sigma, device="cpu")
    else:
        base = build_base_prime_field_numpy(grid_n=args.grid, n_primes=args.n_primes, bump_sigma=args.bump_sigma)

    # Build φ and geodesics at σ=0.5 and off-σ set
    xs = np.linspace(0.12, 0.88, args.seed_grid); ys = np.linspace(0.12, 0.88, args.seed_grid)
    seeds = [(float(x), float(y)) for x in xs for y in ys]
    angles = [k * 2*np.pi/args.angles for k in range(args.angles)]

    def run_sigma(sigma):
        phi = build_phi_from_base(base, sigma, args.eps, args.alpha, args.smooth_px)
        gy, gx = np.gradient(phi)
        errs=[]; tvals=[]
        for (x0,y0) in seeds:
            for th in angles:
                path = geodesic_rk4(x0, y0, th, gx, gy, dt=args.dt, steps=args.steps)
                e, idx = closure_error(path, warmup=160)
                thr = None  # we'll quantile on collected errs later
                win = path[max(0, idx-25):idx+1, 1]
                tvals.append(float(np.mean(win)*args.tmax))
                errs.append(e)
        errs = np.array(errs); tvals = np.array(tvals)
        qthr = float(np.quantile(errs, args.q_closure))
        mask = errs <= qthr
        return {"t": tvals[mask], "errs": errs, "qthr": qthr, "count": int(mask.sum())}

    res_on = run_sigma(0.5)
    res_off = {float(s): run_sigma(s) for s in args.off_sigmas}

    # Predictions and RMSE (t>100)
    zeros = KNOWN_T[KNOWN_T > 100.0]
    pred_on = predict_peaks(res_on["t"], args.tmax, t_min=100.0)
    rmse_on = nearest_zero_rmse(pred_on, zeros)

    # Plot σ=0.5 histogram with zeros and predicted peaks
    plt.figure(figsize=(9,4))
    if len(res_on["t"])>0:
        plt.hist(res_on["t"], bins=80, alpha=0.85)
    for z in zeros: plt.axvline(float(z), linestyle=":", alpha=0.5)
    for c in pred_on: plt.axvline(float(c), linestyle="--", alpha=0.6)
    plt.xlabel("t at σ=0.5"); plt.ylabel("count")
    plt.title("Phase 9: σ=0.5 t‑density with zeros (>100) and predicted peaks")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "phase9_sigma05_t_hist_peaks_gt100.png"), dpi=220)
    plt.close()

    # Off-line comparison: choose nearest off σ (e.g., 0.45 & 0.55) and run permutation null
    perm_results = {}
    for s, res in res_off.items():
        pr = permutation_null_rmse(res_on["t"], res["t"], zeros, R=400, T_MAX=args.tmax, t_min=100.0)
        perm_results[float(s)] = pr

    # Save JSON summary
    summary = {
        "rmse_sigma_05": rmse_on,
        "predicted_peaks_sigma_05": pred_on.tolist(),
        "count_sigma_05": res_on["count"],
        "off_sigma_perm": perm_results
    }
    with open(os.path.join(args.outdir, "phase9_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # Also save raw t-values for on/off
    with open(os.path.join(args.outdir, "phase9_tvals_sigma05.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t_sigma05"]); 
        for v in res_on["t"]: w.writerow([float(v)])
    with open(os.path.join(args.outdir, "phase9_tvals_off.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["sigma","t"]); 
        for s, res in res_off.items():
            for v in res["t"]: w.writerow([float(s), float(v)])

    print("RMSE@σ=0.5 (t>100):", rmse_on)
    print("Permutation null (on - off RMSE) p-values:", perm_results)

if __name__ == "__main__":
    main()
