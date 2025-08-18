#!/usr/bin/env python3
"""Compute RMSE and permutation p-value for predicted peaks against known zeros.

Usage:
    python scripts/validate_metrics.py --pred predicted_peaks.txt --zeros zeros.csv [--perms 1000]

The zeros CSV should have a column containing zero ordinates; the first column is used.
"""

import argparse
import numpy as np


def load_peaks(path: str) -> np.ndarray:
    return np.loadtxt(path, dtype=float)


def load_zeros(path: str) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", skiprows=1, dtype=float)


def compute_rmse(peaks: np.ndarray, zeros: np.ndarray) -> float:
    distances = np.abs(zeros[:, None] - peaks[None, :])
    min_dists = distances.min(axis=0)
    return float(np.sqrt(np.mean(min_dists ** 2)))


def permutation_p(peaks: np.ndarray, zeros: np.ndarray, n_perm: int) -> float:
    obs = compute_rmse(peaks, zeros)
    zmin, zmax = zeros.min(), zeros.max()
    count = 0
    for _ in range(n_perm):
        rand_peaks = np.random.uniform(zmin, zmax, size=len(peaks))
        if compute_rmse(rand_peaks, zeros) <= obs:
            count += 1
    return (count + 1) / (n_perm + 1), obs


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate peak predictions")
    parser.add_argument("--pred", required=True, help="Path to predicted_peaks.txt")
    parser.add_argument("--zeros", required=True, help="Path to zeros CSV")
    parser.add_argument("--perms", type=int, default=1000, help="Number of permutations")
    args = parser.parse_args()

    peaks = load_peaks(args.pred)
    zeros = load_zeros(args.zeros)
    p_val, rmse = permutation_p(peaks, zeros, args.perms)
    print(f"RMSE: {rmse:.6g}")
    print(f"Permutation p-value: {p_val:.6g}")


if __name__ == "__main__":
    main()
