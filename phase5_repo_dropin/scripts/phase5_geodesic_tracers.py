#!/usr/bin/env python3
# (content identical to previous successful definition, shortened for brevity in this comment)
# Please see earlier message for full comments; functional code below.
import argparse, csv
import numpy as np
import matplotlib.pyplot as plt
from math import log

def sieve_primes(n: int):
    sieve = np.ones(n+1, dtype=bool); sieve[:2] = False
    for p in range(2, int(n**0.5)+1):
        if sieve[p]: sieve[p*p:n+1:p] = False
    return [int(i) for i in range(n+1) if sieve[i]]

def bilinear_sample(field: np.ndarray, x: float, y: float) -> float:
    h, w = field.shape; eps = 1e-9
    x = max(0.0, min(1.0 - eps, x)); y = max(0.0, min(1.0 - eps, y))
    gx = x * (w - 1); gy = y * (h - 1)
    x0 = int(np.floor(gx)); y0 = int(np.floor(gy))
    x1 = min(x0 + 1, w - 1); y1 = min(y0 + 1, h - 1)
    tx = gx - x0; ty = gy - y0
    v00 = field[y0, x0]; v10 = field[y0, x1]; v01 = field[y1, x0]; v11 = field[y1, x1]
    return (1-ty)*((1-tx)*v00 + tx*v10) + ty*((1-tx)*v01 + tx*v11)

def bilinear_vec(gx: np.ndarray, gy: np.ndarray, x: float, y: float):
    return bilinear_sample(gx, x, y), bilinear_sample(gy, x, y)

def build_pia_field(grid_n=260, n_primes=160, bump_sigma=0.035):
    field = np.zeros((grid_n, grid_n), dtype=np.float64)
    primes = sieve_primes(900)[:n_primes]
    pmax = primes[-1]; lmax = log(pmax)
    ys, xs = np.mgrid[0:grid_n, 0:grid_n]
    X = xs / (grid_n - 1); Y = ys / (grid_n - 1)
    inv2s2 = 1.0 / (2.0 * bump_sigma * bump_sigma)
    for p in primes:
        xp = p / pmax; yp = log(p) / lmax
        w = np.exp(-((X - xp)**2 + (Y - yp)**2) * inv2s2) / p
        field += w
    return np.pi + field

def smooth_gaussian(field: np.ndarray, sigma_px=1.8):
    radius = int(3*sigma_px); xs = np.arange(-radius, radius+1)
    ker = np.exp(-0.5*(xs/sigma_px)**2); ker /= ker.sum()
    tmp = np.apply_along_axis(lambda m: np.convolve(m, ker, mode="same"), axis=1, arr=field)
    return np.apply_along_axis(lambda m: np.convolve(m, ker, mode="same"), axis=0, arr=tmp)

def integrate_flow(x0, y0, gx, gy, steps=900, dt=0.008):
    pos = np.array([x0, y0], dtype=np.float64); path = [pos.copy()]
    def F(p):
        px = max(0.0, min(1.0, float(p[0]))); py = max(0.0, min(1.0, float(p[1])))
        fx, fy = bilinear_vec(gx, gy, px, py); v = np.array([-fx, -fy])
        n = np.linalg.norm(v) + 1e-9; return v / n
    for _ in range(steps):
        k1 = F(pos); k2 = F(pos + 0.5*dt*k1); k3 = F(pos + 0.5*dt*k2); k4 = F(pos + dt*k3)
        pos = pos + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        pos[0] = min(max(pos[0], 0.0), 1.0); pos[1] = min(max(pos[1], 0.0), 1.0)
        path.append(pos.copy())
    return np.array(path)

def closure_error(path: np.ndarray, min_period=160):
    start = path[0]; seg = path[min_period:]
    if len(seg) == 0: return 1e9, min_period
    d = np.linalg.norm(seg - start, axis=1); i = int(np.argmin(d))
    return float(d[i]), min_period + i

def main():
    import argparse, os
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--grid", type=int, default=260)
    ap.add_argument("--primes", type=int, default=160)
    ap.add_argument("--steps", type=int, default=900)
    ap.add_argument("--dt", type=float, default=0.008)
    ap.add_argument("--closure_tol", type=float, default=0.035)
    ap.add_argument("--bump_sigma", type=float, default=0.035)
    ap.add_argument("--smooth_px", type=float, default=2.0)
    ap.add_argument("--tmax", type=float, default=50.0, help="map y∈[0,1] → t∈[0,T_MAX]")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    known_t = np.array([14.13472514, 21.02203964, 25.01085758, 30.42487613, 32.93506159,
                        37.58617816, 40.91871901, 43.32707328, 48.00515088, 49.77383248])

    pia = build_pia_field(grid_n=args.grid, n_primes=args.primes, bump_sigma=args.bump_sigma)
    pia = smooth_gaussian(pia, sigma_px=args.smooth_px)
    phi = (pia - np.pi); phi /= (np.max(np.abs(phi)) + 1e-9)
    gx, gy = np.gradient(phi)

    seeds = [(float(x0), float(y0)) for x0 in np.linspace(0.1,0.9,14) for y0 in np.linspace(0.1,0.9,14)]
    closures, paths_subset = [], []
    for (x0, y0) in seeds:
        path = integrate_flow(x0, y0, gx, gy, steps=args.steps, dt=args.dt)
        err, idx = closure_error(path, min_period=160)
        if err < args.closure_tol:
            win = path[max(0, idx-30):idx+1, 1]; t_est = float(np.mean(win) * args.tmax)
            closures.append((x0, y0, err, t_est))
            if len(paths_subset) < 60: paths_subset.append(path)

    import csv
    with open(os.path.join(args.outdir, "closure_matches.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["x0","y0","closure_error","t_est","nearest_zero_dist"])
        for x0, y0, err, t_est in closures:
            nz = float(np.min(np.abs(known_t - t_est))); w.writerow([x0, y0, err, t_est, nz])

    plt.figure(figsize=(6,6)); plt.imshow(phi, origin='lower', extent=[0,1,0,1], aspect='equal', cmap="viridis")
    for path in paths_subset: plt.plot(path[:,0], path[:,1], linewidth=0.7, alpha=0.9)
    plt.title("πₐ field (φ) with geodesic flow trajectories")
    plt.xlabel("x (normalized)"); plt.ylabel("y (normalized)  ~  t / T_MAX"); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "pia_field_trajectories.png"), dpi=220); plt.close()

    t_ests = np.array([t for _,_,_,t in closures]) if closures else np.array([])
    plt.figure(figsize=(9,4))
    if len(t_ests) > 0: plt.hist(t_ests, bins=30, alpha=0.8, label="closure t_est")
    for z in known_t: plt.axvline(z, linestyle="--", alpha=0.5, label="known zeros" if z==known_t[0] else None)
    plt.xlabel("t (imaginary part estimate from closures)"); plt.ylabel("count")
    plt.title("Closure-based zero estimates vs known zeros"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "closure_t_hist.png"), dpi=220); plt.close()

if __name__ == "__main__":
    main()
