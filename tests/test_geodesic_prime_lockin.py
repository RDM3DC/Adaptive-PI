import math
import importlib.util
import types
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
pkg = types.ModuleType("pi_a")
pkg.__path__ = [str(_SRC / "pi_a")]
sys.modules.setdefault("pi_a", pkg)

spec_m = importlib.util.spec_from_file_location("pi_a.models", _SRC / "pi_a" / "models.py")
models = importlib.util.module_from_spec(spec_m)
spec_m.loader.exec_module(models)
sys.modules.setdefault("pi_a.models", models)

spec_g = importlib.util.spec_from_file_location("pi_a.geodesics", _SRC / "pi_a" / "geodesics.py")
geodesics = importlib.util.module_from_spec(spec_g)
spec_g.loader.exec_module(geodesics)
sys.modules.setdefault("pi_a.geodesics", geodesics)

geodesic_flow = geodesics.geodesic_flow
ConstantCurvature = models.ConstantCurvature


def test_prime_lockin_behaviour():
    steps = 200
    step = 0.05
    # curvature chosen so that without prime forcing the path closes
    K0 = ConstantCurvature(2 * math.pi / (steps * step))

    path = geodesic_flow((0.0, 0.0), 0.0, K0, steps=steps, step=step,
                         prime_weight=1.0, re_s=0.5)
    x_end, y_end = path[-1]
    d1 = math.hypot(x_end, y_end)
    assert d1 < 0.2

    path2 = geodesic_flow((0.0, 0.0), 0.0, K0, steps=steps, step=step,
                          prime_weight=1.0, re_s=0.6)
    x2, y2 = path2[-1]
    d2 = math.hypot(x2, y2)
    assert d2 > d1 + 0.2
