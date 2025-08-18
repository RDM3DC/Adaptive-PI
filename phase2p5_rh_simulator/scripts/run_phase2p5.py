#!/usr/bin/env python3
"""
Phase II.5 — Adaptive π (π_a) RH Geodesic Simulator
- Weighted π_a by 1/(p * |σ - 0.5| + ε)
- Gaussian smoothing (no SciPy required)
- Conformal-metric geodesics via Christoffel form (2D, e^{2φ} δ_ij)
- Outputs: closure vs σ, φ-field snapshot, geodesic overlay, closure heatmap
"""

import numpy as np
import matplotlib.pyplot as plt
from math import log
import json
from typing import List, Tuple

# ------------------ Prime utilities ------------------
def sieve_primes(n: int) -> List[int]:
    sieve = np.ones(n+1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(n**0.5)+1):
        if sieve[p]:
            sieve[p*p:n+1:p] = False
    return [int(i) for i in range(n+1) if sieve[i]]

# ------------------ Image helpers ------------------
def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()

# ------------------ Field sampling ------------------
def bilinear_sample(field: np.ndarray, x: float, y: float) -> float:
    h, w = field.shape
    gx = x * (w - 1)
    gy = y * (h - 1)
    x0 = int(np.floor(gx)); y0 = int(np.floor(gy))
    x1 = min(x0 + 1, w - 1); y1 = min(y0 + 1, h - 1)
    tx = gx - x0; ty = gy - y0
    v00 = field[y0, x0]; v10 = field[y0, x1]; v01 = field[y1, x0]; v11 = field[y1, x1]
    return (1-ty)*((1-tx)*v00 + tx*v10) + ty*((1-tx)*v01 + tx*v11)

# ------------------ π_a field ------------------
def build_pi_a_field(sigma: float, grid_n: int = 300, n_primes: int = 150, bump_sigma: float = 0.035, eps: float = 1e-3):
    field = np.zeros((grid_n, grid_n), dtype=np.float64)
    primes = sieve_primes(1000)[:n_primes]
    p_max = primes[-1]; log_p_max = log(p_max)
    bias = 1.0 / (abs(sigma - 0.5) + eps)

    ys, xs = np.mgrid[0:grid_n, 0:grid_n]
    X = xs / (grid_n - 1); Y = ys / (grid_n - 1)
    inv2s2 = 1.0 / (2.0 * bump_sigma * bump_sigma)

    for p in primes:
        xp = p / p_max; yp = log(p) / log_p_max
        dx2 = (X - xp) ** 2; dy2 = (Y - yp) ** 2
        w = np.exp(-(dx2 + dy2) * inv2s2)
        field += (w / p) * bias

    pi_a = np.pi + field
    return pi_a

def smooth_field(field: np.ndarray, sigma_px: float = 2.0) -> np.ndarray:
    radius = int(3 * sigma_px)
    xs = np.arange(-radius, radius+1)
    kernel = np.exp(-0.5 * (xs / sigma_px) ** 2); kernel /= kernel.sum()
    tmp = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=1, arr=field)
    smoothed = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=0, arr=tmp)
    return smoothed

def grad_field(field: np.ndarray):
    gy, gx = np.gradient(field)
    return gx, gy

# ------------------ Geodesic integrator ------------------
def geodesic_rk4(x0, y0, theta, phi, gx_field, gy_field, dt=0.004, steps=1200):
    v = np.array([np.cos(theta), np.sin(theta)], dtype=np.float64)
    pos = np.array([x0, y0], dtype=np.float64)
    path = [pos.copy()]

    def grad_phi_at(px, py):
        return np.array([bilinear_sample(gx_field, px, py), bilinear_sample(gy_field, px, py)], dtype=np.float64)

    for _ in range(steps):
        def acc(p, v):
            g = grad_phi_at(p[0], p[1])
            vv = np.dot(v, v)
            return -2.0 * np.dot(g, v) * v + g * vv

        p1 = pos; v1 = v; a1 = acc(p1, v1)
        p2 = pos + 0.5*dt*v1; v2 = v + 0.5*dt*a1; a2 = acc(p2, v2)
        p3 = pos + 0.5*dt*v2; v3 = v + 0.5*dt*a2; a3 = acc(p3, v3)
        p4 = pos + dt*v3;     v4 = v + dt*a3;     a4 = acc(p4, v4)

        pos = pos + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
        v   = v   + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

        if pos[0] < 0.0 or pos[0] > 1.0:
            v[0] *= -1.0; pos[0] = np.clip(pos[0], 0.0, 1.0)
        if pos[1] < 0.0 or pos[1] > 1.0:
            v[1] *= -1.0; pos[1] = np.clip(pos[1], 0.0, 1.0)
        path.append(pos.copy())
    return np.array(path)

def closure_error(path, min_period=200):
    start = path[0]
    seg = path[min_period:]
    dists = np.linalg.norm(seg - start, axis=1)
    return float(np.min(dists))

# ------------------ Experiments ------------------
def run_sigma_sweep(outdir="outputs", seed_points=None, angles=None, sigmas=None):
    os.makedirs(outdir, exist_ok=True)
    if seed_points is None:
        seed_points = [(0.2,0.2), (0.8,0.2), (0.2,0.8), (0.8,0.8), (0.5,0.5), (0.3,0.6)]
    if angles is None:
        angles = [0.0, np.pi/3, 2*np.pi/3, np.pi, 4*np.pi/3]
    if sigmas is None:
        sigmas = np.linspace(0.15, 0.85, 11)

    results = {}
    for s in sigmas:
        pi_a = build_pi_a_field(float(s))
        pi_a = smooth_field(pi_a, sigma_px=2.0)
        centered = pi_a - np.pi
        centered /= (np.max(np.abs(centered)) + 1e-9)
        phi = 0.4 * centered
        gx, gy = grad_field(phi)

        errs = []
        for (x0,y0) in seed_points:
            for th in angles:
                path = geodesic_rk4(x0, y0, th, phi, gx, gy, dt=0.004, steps=1100)
                errs.append(closure_error(path, min_period=220))

        results[float(s)] = float(np.mean(errs))

    # Save results + figure
    with open(os.path.join(outdir, "closure_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    xs = list(results.keys()); ys = list(results.values())
    plt.figure(figsize=(7,4))
    plt.plot(xs, ys, marker='o')
    plt.axvline(0.5, linestyle='--')
    plt.xlabel("σ = Re(s)"); plt.ylabel("Mean closure error")
    plt.title("Phase II.5: Closure vs Re(s)")
    plt.grid(True, linestyle=":")
    savefig(os.path.join(outdir, "closure_vs_sigma.png"))

    # Field + geodesics at σ=0.5
    s_focus = 0.5
    pi_a = build_pi_a_field(s_focus)
    pi_a = smooth_field(pi_a, sigma_px=2.0)
    centered = pi_a - np.pi; centered /= (np.max(np.abs(centered)) + 1e-9)
    phi = 0.4 * centered
    gx, gy = grad_field(phi)

    plt.figure(figsize=(6,6))
    plt.imshow(phi, origin='lower', extent=[0,1,0,1], aspect='equal')
    plt.title("φ field (σ = 0.5)")
    plt.xlabel("x"); plt.ylabel("y")
    savefig(os.path.join(outdir, "field_sigma05.png"))

    plt.figure(figsize=(6,6))
    plt.imshow(phi, origin='lower', extent=[0,1,0,1], aspect='equal')
    for (x0,y0) in seed_points:
        for th in angles:
            path = geodesic_rk4(x0, y0, th, phi, gx, gy, dt=0.004, steps=900)
            plt.plot(path[:,0], path[:,1], linewidth=0.9)
    plt.title("Geodesics over φ (σ = 0.5)")
    plt.xlabel("x"); plt.ylabel("y")
    savefig(os.path.join(outdir, "geodesics_sigma05.png"))

    # Closure heatmap at σ=0.5
    grid = np.linspace(0.1, 0.9, 22)
    HH = np.zeros((len(grid), len(grid)))
    for i, y0 in enumerate(grid):
        for j, x0 in enumerate(grid):
            path = geodesic_rk4(float(x0), float(y0), 0.0, phi, gx, gy, dt=0.004, steps=650)
            HH[i,j] = closure_error(path, min_period=180)

    plt.figure(figsize=(6,5.5))
    plt.imshow(HH, origin='lower', extent=[0.1,0.9,0.1,0.9], aspect='auto')
    plt.colorbar(label="Closure error")
    plt.title("Closure error heatmap (σ = 0.5)")
    plt.xlabel("x0"); plt.ylabel("y0")
    savefig(os.path.join(outdir, "closure_heatmap_sigma05.png"))

if __name__ == "__main__":
    import os, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="outputs")
    parser.add_argument("--grid", type=int, default=300, help="internal grid size (default 300)")
    parser.add_argument("--primes", type=int, default=150, help="number of primes (default 150)")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    run_sigma_sweep(outdir=args.outdir)
