<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepCurveExt: curves

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Curve operations](#curve-operations)
- [Coincidence of projected curves](#coincidence-of-projected-curves)
- [Stock counterpart](#stock-counterpart)

## Notes by topic

### Curve operations

Offsets of planar curves with three corner treatments (arc, chord, extension to the crossing) are made on a
copy moved to the origin and scaled, so OCCT's absolute 1e-7 lengths do not change the result with the curve's size.
The outline of a body along a direction uses the edges and the turning curve `n.v = 0` of curved faces, with the
places where visibility can change found exactly. Tangent circles are found by Newton in a frame of the curves' own
extent. Also: projection, rebuild with a prescribed degree and count, refinement, and regions of coplanar curves.
A spiral is the helix (angle `2 pi t`, radius `R + p t`, height `h t / turns`).

### Coincidence of projected curves

Five common points and equal lengths do not imply equal curves. Plane projection removes an edge only after
establishing a common normalized spline basis (degree elevation, knot insertion, optionally reversal) with identical
positive weights. If every corresponding pole differs by at most the modelling tolerance, the rational
convex-combination formula bounds the displacement at every parameter by that tolerance, in both directions since the
common interval covers both edges, with outward-rounded pole-distance arithmetic [PMC09 §1.4.2] (derived). Edges whose
correspondence cannot be established are both retained. Degree elevation and knot insertion are OCCT's floating-point
implementations; the enclosure bounds the resulting common-basis representations, with the usual rounding limits.

### Stock counterpart

Offsets: `BRepOffsetAPI_MakeOffset`. Outline: `HLRBRep` (its documented limits on tangent and edge-on faces
are what the exact event search avoids). The rest has no direct stock counterpart.
