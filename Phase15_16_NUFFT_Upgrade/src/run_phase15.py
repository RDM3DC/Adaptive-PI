
import argparse, os, json, numpy as np
from src.prime_gen import first_n_primes
from src.nodes_and_weights import log_nodes_from_primes, weights_from_primes, make_t_grid
from src.nufft_eval import eval_many_t
from src.metrics import rmse_to_zeros

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_primes", type=int, default=1_000_000)
    ap.add_argument("--tmax", type=float, default=21000.0)
    ap.add_argument("--tstep", type=float, default=0.25)
    ap.add_argument("--method", type=str, default="finufft", choices=["finufft","torch","naive"])
    ap.add_argument("--zeros_csv", type=str, default="")
    ap.add_argument("--tmin_eval", type=float, default=20000.0)
    ap.add_argument("--outdir", type=str, default="outputs_phase15")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"[phase15] generating {args.n_primes} primes…")
    primes, bound = first_n_primes(args.n_primes)
    primes = np.array(primes, dtype=np.uint64)

    print(f"[phase15] building nodes/weights…")
    x = log_nodes_from_primes(primes)
    w = weights_from_primes(primes)

    print(f"[phase15] building t-grid up to {args.tmax}…")
    t = make_t_grid(0.0, args.tmax, args.tstep)

    print(f"[phase15] NUFFT method={args.method} …")
    H = eval_many_t(x, w, t, method=args.method)

    # Simple trough detection via second-derivative sign changes
    # (Use your more sophisticated peak picker in production.)
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
        "n_primes": int(args.n_primes),
        "tmax": float(args.tmax),
        "tstep": float(args.tstep),
        "method": args.method,
        "rmse_tmin": float(args.tmin_eval),
        "rmse": None if rmse is None else float(rmse),
    }
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print("[phase15] done; wrote:", args.outdir)

if __name__ == "__main__":
    main()
