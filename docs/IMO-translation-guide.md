# Translating IMO Geometry to \u03c0\u2090

This guide shows how to port classic Euclidean arguments to \u03c0\u2090.

1) **Replace** any use of \u201c\u03c0\u201d from straight/half-turn closures with **\u03c0\u2090(region)** = \u03c0 + \u03a6_region.
2) **Angle-chasing:** carry a flux budget. When two inscribed angles are compared, add the sector-flux difference.
3) **Right angles:** remain right in the flat limit; away from flatness, use local orthogonality modulo O(\u03a6) corrections.
4) **Ceva/Menelaus:** use trig forms and attach exp(\u03b5\u00b7\u039e) correction, which vanishes when K\u22610.
5) **Inversion/circle power:** approximate by \u03c0\u2090-circles (geodesic circles). To first order, radii/centers perturb; equal-power relations include small flux terms but preserve structure.

Worked examples are in `examples/`.
