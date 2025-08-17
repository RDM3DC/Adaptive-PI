from __future__ import annotations
from typing import Tuple, List
import math

from .models import CurvatureField

Point = Tuple[float, float]


def geodesic_circle(O: Point, r: float, K: CurvatureField, n: int = 128) -> List[Point]:
    """Sample a 'geodesic circle' using Euclidean proxy; for small curvature treat as Euclidean circle.
    Future refinement: adjust radial distance based on average curvature along rays.
    """
    pts = []
    for i in range(n):
        th = 2 * math.pi * i / n
        pts.append((O[0] + r * math.cos(th), O[1] + r * math.sin(th)))
    return pts


def _is_prime(n: int) -> bool:
    """Return ``True`` if *n* is a prime number."""
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True


def geodesic_flow(
    start: Point,
    angle: float,
    K: CurvatureField,
    steps: int = 100,
    step: float = 0.1,
    prime_weight: float = 0.0,
    re_s: float = 0.5,
) -> List[Point]:
    """Integrate a geodesic under a πₐ curvature field.

    Parameters
    ----------
    start:
        Starting point ``(x, y)``.
    angle:
        Initial heading in radians.
    K:
        Curvature field used to bend the trajectory.
    steps:
        Number of integration steps.
    step:
        Step size for each integration increment.
    prime_weight:
        Extra curvature applied when the step index is prime. The amount is
        scaled by ``(re_s - 1/2)`` so that the forcing vanishes on the
        critical line ``Re(s) = 1/2``.
    re_s:
        Real part of the zeta argument modelling distance from the critical
        line.

    Returns
    -------
    list[Point]
        Sampled points along the flow, including the starting point.
    """

    pts = [start]
    x, y = start
    th = angle
    for i in range(1, steps + 1):
        k = K(x, y)
        if prime_weight and _is_prime(i):
            k += prime_weight * (re_s - 0.5)
        th += k * step
        x += step * math.cos(th)
        y += step * math.sin(th)
        pts.append((x, y))
    return pts
