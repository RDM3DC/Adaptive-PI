
import numpy as np, math
from src.nodes_and_weights import make_t_grid
# small primes for quick run
def k_primes(k):
    limit = 200000
    sieve = np.ones(limit+1, dtype=bool); sieve[:2]=False
    for p in range(2,int(limit**0.5)+1):
        if sieve[p]: sieve[p*p:limit+1:p]=False
    ps = np.nonzero(sieve)[0][:k]
    return ps
ps = k_primes(10000)
x = np.log(ps.astype(float))
w = ps.astype(float)**(-0.5)
t = make_t_grid(0.0, 200.0, 0.25)
H = np.zeros_like(t)
for j in range(len(ps)):
    H += np.cos(t * x[j]) * w[j]
np.save("demo/H.npy", H)
np.save("demo/t.npy", t)
print("Demo wrote demo/H.npy and demo/t.npy")
