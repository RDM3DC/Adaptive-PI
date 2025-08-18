
import argparse, os, glob, json
import numpy as np
from typing import List
from src.nufft_eval import eval_many_t
from src.nodes_and_weights import make_t_grid
from src.metrics import rmse_to_zeros

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmax", type=float, required=True)
    ap.add_argument("--tstep", type=float, default=0.20)
    ap.add_argument("--method", type=str, default="finufft")
    ap.add_argument("--zeros_csv", type=str, default="")
    ap.add_argument("--tmin_eval", type=float, default=30000.0)
    ap.add_argument("--shard_glob", type=str, required=True,
                    help="glob for shards, each .npy contains log(p) as float64")
    ap.add_argument("--outdir", type=str, default="outputs_phase16")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    t = make_t_grid(0.0, args.tmax, args.tstep)
    H = np.zeros_like(t, dtype=np.float64)

    shards = sorted(glob.glob(args.shard_glob))
    if not shards:
        raise SystemExit(f"No shards found for pattern: {args.shard_glob}")

    for i, shard_path in enumerate(shards, 1):
        x = np.load(shard_path)  # log(p) values
        w = np.exp(-0.5 * x)     # since p = exp(x), 1/sqrt(p) = e^{-x/2}
        print(f"[phase16] shard {i}/{len(shards)}: {os.path.basename(shard_path)} (N={len(x)})")
        H += eval_many_t(x, w.astype(np.complex128), t, method=args.method)
        # optional checkpoint
        np.save(os.path.join(args.outdir, f"H_partial_{i:04d}.npy"), H)

    # crude troughs (replace with your peak detector if desired)
    d1 = np.gradient(H, args.tstep)
    d2 = np.gradient(d1, args.tstep)
    trough_idx = np.where((d1[1:-1]*d1[2:] <= 0) & (d2[1:-1] > 0))[0] + 1
    pred_t = t[trough_idx]

    rmse = None
    if args.zeros_csv and os.path.exists(args.zeros_csv):
        zeros = np.loadtxt(args.zeros_csv, dtype=float)
        rmse = rmse_to_zeros(pred_t, zeros, t_min=args.tmin_eval)

    np.save(os.path.join(args.outdir, "t_grid.npy"), t)
    np.save(os.path.join(args.outdir, "H.npy"), H)
    np.savetxt(os.path.join(args.outdir, "predicted_peaks.txt"), pred_t, fmt="%.6f")
    meta = {
        "tmax": float(args.tmax),
        "tstep": float(args.tstep),
        "method": args.method,
        "tmin_eval": float(args.tmin_eval),
        "rmse": None if rmse is None else float(rmse),
        "num_shards": len(shards),
        "shard_glob": args.shard_glob,
    }
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print("[phase16] done; wrote:", args.outdir)

if __name__ == "__main__":
    main()
