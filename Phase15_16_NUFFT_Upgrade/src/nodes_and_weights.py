
import numpy as np
from typing import Tuple, Optional

def weights_from_primes(primes: np.ndarray, weight: float = 0.5) -> np.ndarray:
    return (primes.astype(np.float64))**(-weight)

def log_nodes_from_primes(primes: np.ndarray) -> np.ndarray:
    return np.log(primes.astype(np.float64))

def make_t_grid(t_min: float, t_max: float, t_step: float) -> np.ndarray:
    return np.arange(t_min, t_max + 1e-12, t_step)

def shard_iterable(arr: np.ndarray, shard_size: int):
    n = len(arr)
    for i in range(0, n, shard_size):
        yield arr[i:i+shard_size]
