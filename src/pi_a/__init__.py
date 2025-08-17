try:
    from .core import (
        triangle_angle_sum_pia,
        polygon_angle_sum_pia,
        sector_flux,
        is_pia_cyclic,
    )
except Exception:  # pragma: no cover - optional deps like numpy may be absent
    triangle_angle_sum_pia = polygon_angle_sum_pia = sector_flux = is_pia_cyclic = None

from .models import ConstantCurvature, GaussianBump
