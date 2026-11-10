<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepBooleanExt: cell-based Booleans

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Booleans](#booleans)
- [Unproven choices](#unproven-choices)

## Notes by topic

### Booleans

Cell selection from one General Fuse with a ladder of fallbacks, each rung accepted only if checked; a sheet tool
divides rather than subtracts.

For a linear extrusion `S(u,v)=C(u)+vD`, signed volume density is
`g(u)=(C(u)-O).(C'(u) cross D)/3`, independent of `v`. Green's theorem reduces
the face contribution to `-integral (v-v0)g(u) du` around its oriented parameter
wires. Holes keep their winding and reversed faces reverse the result. The
construction permits boundaries that fold over `u` and unwrapped periodic
charts. Stored pcurve endpoints can differ within edge tolerance; straight
parameter segments close every consecutive pair before integration. Without
these segments the omitted flux depends on the primitive's origin.

For a polynomial carrier of degree `p`, the leading determinant coefficient
vanishes, so `degree(g)<=2p-2`. Along a degree-`q` polynomial trim, the contour
integrand has degree at most `2pq-1`: `pq` Gauss-Legendre nodes therefore
integrate it exactly in real arithmetic. Integration splits at trim knots and
pullbacks of carrier knots. Closing segments require `p` nodes. Rectangles
retain the existing eight-point rule when `p<=8`, where it is already exact;
higher degrees use `p` nodes per carrier span. Degree-exact Legendre projections, converted to finite Chebyshev expansions,
cache the flux and trim jets; Clenshaw evaluation avoids
repeated spline evaluation. Two-double arithmetic preserves the stored binary64
geometry through evaluation and accumulation before the final binary64 result.

Rational carriers or trims retain the existing adaptive property integrator,
including its Gauss-Kronrod rule, tolerance, stopping and convergence policy.
Area, moments, the coarse material bracket, Boolean tolerance and allowance
formulas are unchanged. Only the integrated flux in the quadrature allowance
can change numerically. Polynomial exactness is subject to floating-point
rounding and computed knot intersections; neither this result nor the adaptive
integration provides an interval-certified error bound.

| Analytic regression | Carrier degree | Trim degree | Gauss nodes | Absolute error limit |
|---|---:|---:|---:|---:|
| Curved boundary, both orientations and tolerance-sized endpoint gap | 2 | 2 | 4 | 3e-16 |
| Rectangle with degree-16 flux | 9 | 1 | 9 | 2e-17 |
| Periodic square carrier, three charts with and without endpoint gaps | 1 | 1 | 1 per span | exact binary64 result |

The dimensional reduction follows Gunderman, Weiss and Evans,
*High-accuracy mesh-free quadrature for trimmed parametric surfaces and volumes*
(2021), sections 2.1.2 and 2.2. The extrusion degree bounds above are derived
from the determinant. Two-double operations follow Hida, Li and Bailey,
*Algorithms for Quad-Double Precision Floating Point Arithmetic* (2000), section 2.

### Unproven choices

Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean
fallbacks, the sampling counts of the offset and skin validity checks.
