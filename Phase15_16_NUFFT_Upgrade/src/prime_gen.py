
# NOTE: For serious runs, install primesieve (pip install primesieve) and use its fast generator.
# Fallback segmented sieve is provided for portability.

from typing import Tuple, Iterator, List
import math

try:
    import primesieve as _ps
except Exception:
    _ps = None

def simple_sieve(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            start = p*p
            step = p
            sieve[start:limit+1:step] = b"\x00" * (((limit - start)//step) + 1)
    return [i for i in range(limit + 1) if sieve[i]]

def segmented_sieve(limit: int, segment_size: int = 1<<21):
    small = simple_sieve(int(limit**0.5) + 1)
    for p in small:
        yield p
    low = max(2, int(limit**0.5) + 1)
    while low <= limit:
        high = min(low + segment_size - 1, limit)
        size = high - low + 1
        mark = bytearray(b"\x01") * size
        for p in small:
            start = (-(low % p)) % p
            start = start if (low + start) >= p*p else (p*p - low)
            if start < 0: start = 0
            for j in range(start, size, p):
                mark[j] = 0
        for i in range(size):
            if mark[i]:
                yield low + i
        low = high + 1

def first_n_primes(n: int) -> Tuple[list, int]:
    if _ps is not None:
        # fast path: primesieve count + nth_prime
        # we stream until n reached
        out = []
        it = _ps.Iterator()
        p = it.next_prime()
        while p and len(out) < n:
            out.append(p)
            p = it.next_prime()
        return out[:n], out[-1]
    # fallback approximate limit
    if n < 6:
        limit = 20
    else:
        x = n * (math.log(n) + math.log(math.log(n)))
        limit = int(x + 10*n)
    out = []
    while True:
        out.clear()
        for p in segmented_sieve(limit):
            out.append(p)
            if len(out) >= n:
                return out[:n], limit
        limit = int(limit * 1.6)
