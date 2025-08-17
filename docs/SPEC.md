# \u03c0\u2090 Geometry \u2014 Reference Spec (v0.1)

## 1. Principle
Let K(x,y) denote Gaussian curvature (or an effective curvature density) over a planar domain D. For any geodesic polygon P with interior D_P,

**Triangle angle sum (first order):**
S_\u0394 := A + B + C = \u03c0 + \u03a6_\u0394,   where  \u03a6_\u0394 := \u222c\u222c_{D_\u0394} K dA.

**m-gon angle sum (first order):**
S_m = (m-2)\u03c0 + \u03a6_P, with \u03a6_P := \u222c\u222c_{D_P} K dA.

**Cyclicity (adaptive inscribed-angle):**
In Euclidean space, A,B,C,X are concyclic iff \u2220ABX = \u2220ACX. In \u03c0\u2090, we model the inscribed angle as picking up a sector flux correction:

\u2220ABX - \u2220ACX \u2248 \u00bd(\u03a6_{ABX}^\u25D6 \u2212 \u03a6_{ACX}^\u25D6).

A,B,C,X are **\u03c0\u2090-cyclic** if the difference above vanishes to first order (sectors balance), recovering the Euclidean criterion when K\u22610.

**Ceva (trig form, \u03c0\u2090):**
For cevians AD, BE, CF concurrent in triangle ABC,

( sin(\u00c2_D) sin( \u0042\u0302_E ) sin( \u00c7_F ) ) / ( sin(\u00c2'_D) sin( \u0042\u0302'_E ) sin( \u00c7'_F ) ) = 1 \u00b7 exp(\u03b5\u00b7\u039e) ,

where the numerator/denominator use the appropriate interior angles at each split and \u039e is a linear functional of local flux imbalances (vanishes in K\u22610). For small curvature, exp(\u03b5\u00b7\u039e) \u2248 1 + \u03b5\u00b7\u039e.

This spec fixes signs and orientations in code; see `core.py` for definitions.

## 2. Flat limit
All formulas reduce to Euclidean identities when \u03a6 \u2261 0.

## 3. Validity
First-order in |K|\u00b7Area and for slowly varying K across the polygon/sector footprint. For larger curvature or strong gradients, numerical geodesic integrations are required (future work).
