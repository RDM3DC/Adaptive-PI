from __future__ import annotations
import math
from typing import Callable

class CurvatureField:
    def __call__(self, x: float, y: float) -> float:
        raise NotImplementedError

class ConstantCurvature(CurvatureField):
    def __init__(self, K: float):
        self.K = K
    def __call__(self, x: float, y: float) -> float:
        return self.K

class GaussianBump(CurvatureField):
    def __init__(self, K0: float, x0: float, y0: float, sigma: float):
        self.K0, self.x0, self.y0, self.sigma = K0, x0, y0, sigma
    def __call__(self, x: float, y: float) -> float:
        dx, dy = x - self.x0, y - self.y0
        r2 = (dx*dx + dy*dy)
        return self.K0 * math.exp(-r2/(2*self.sigma*self.sigma))
