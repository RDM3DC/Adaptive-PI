
import numpy as np
from typing import Dict

def rmse_to_zeros(pred_t: np.ndarray, zeros: np.ndarray, t_min: float = 0.0) -> float:
    pred = pred_t[pred_t >= t_min]
    if pred.size == 0 or zeros.size == 0: return float("nan")
    d = np.min(np.abs(pred[:,None] - zeros[None,:]), axis=1)
    return float(np.sqrt(np.mean(d**2)))

def permutation_pvalue(on: np.ndarray, off: np.ndarray, zeros: np.ndarray,
                       R: int = 5000, t_min: float = 0.0, seed: int = 0) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    pred_on = on[on>=t_min]; pred_off = off[off>=t_min]
    def _rmse(pred):
        if pred.size == 0 or zeros.size == 0: return np.inf
        d = np.min(np.abs(pred[:,None] - zeros[None,:]), axis=1)
        return float(np.sqrt(np.mean(d**2)))
    rm_on = _rmse(pred_on); rm_off = _rmse(pred_off)
    pool = np.concatenate([pred_on, pred_off])
    n_on = pred_on.size
    diffs = []
    for _ in range(R):
        rng.shuffle(pool)
        a = pool[:n_on]; b = pool[n_on:]
        diffs.append(_rmse(a) - _rmse(b))
    diffs = np.array(diffs)
    obs = rm_on - rm_off
    p = (np.sum(np.abs(diffs) >= abs(obs)) + 1) / (len(diffs) + 1)
    return {"rmse_on": rm_on, "rmse_off": rm_off, "p": float(p), "obs_diff": float(obs)}
