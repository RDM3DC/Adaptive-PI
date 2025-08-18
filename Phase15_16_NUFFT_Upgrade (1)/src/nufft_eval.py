

import numpy as np
from typing import Optional, Dict

# We treat H(t) = Re \sum_j w_j exp(i t x_j)
# where x_j = log p_j (nodes), w_j = 1/sqrt(p_j)
# We support three modes:
#  - FINUFFT type-1 (preferred): pip install finufft (CPU) or cuFINUFFT (GPU)
#  - Torch-NUFFT fallback if available (not exact type-1)
#  - Naive block evaluation for small N (for smoke tests)

def eval_many_t(nodes_x: np.ndarray,
                weights: np.ndarray,
                t_vals: np.ndarray,
                method: str = "finufft",
                iflag: int = +1,
                eps: float = 1e-12,
                block_n: int = 2_000_000) -> np.ndarray:
    nodes_x = nodes_x.astype(np.float64)
    weights = weights.astype(np.complex128)  # complex for exp
    t_vals  = t_vals.astype(np.float64)

    if method == "finufft":
        try:
            import finufft
            # Scale to type-1: nonuniform x -> uniform modes k
            # Map t to integer modes by scaling: k = round(t * a), choose a so that k grid fits desired t grid
            # Simpler: use FINUFFT's type-3 (nonuniform x and nonuniform freq), available in finufft v2.x
            # If type-3 available:
            if hasattr(finufft, "nufft3"):
                fk = finufft.nufft3(type=1, x=nodes_x, y=None, z=None, c=weights,
                                    s=t_vals, t=None, u=None, isign=iflag, eps=eps)
                return np.real(fk)
            # Else fallback to block naive (shouldn't happen often)
        except Exception:
            pass

    if method == "torch":
        try:
            import torch
            X = torch.tensor(nodes_x)[None,:]      # [1, N]
            W = torch.tensor(weights)              # [N]
            out = []
            for tb in torch.tensor(t_vals).split(4096):
                T = tb[:,None]                     # [B,1]
                z = torch.exp(1j * T * X).matmul(W)  # [B,N]@[N]->[B]
                out.append(z.detach().cpu().numpy())
            return np.real(np.concatenate(out, axis=0))
        except Exception:
            pass

    # Naive block evaluation (vectorized per block)
    vals = np.zeros_like(t_vals, dtype=np.complex128)
    for i in range(0, len(nodes_x), block_n):
        x = nodes_x[i:i+block_n]
        w = weights[i:i+block_n]
        # outer: t_vals x x -> exp(i * t * x)
        E = np.exp(1j * t_vals[:,None] * x[None,:])
        vals += (E @ w)
    return np.real(vals)
