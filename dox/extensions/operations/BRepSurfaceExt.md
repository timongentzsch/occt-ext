<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepSurfaceExt: surfaces

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Surface operations](#surface-operations)
- [Alignment input-shape prior](#alignment-input-shape-prior)
- [Collinear alignment: units and rational boundary rows](#collinear-alignment-units-and-rational-boundary-rows)
- [Rectangular sheet extension](#rectangular-sheet-extension)
- [Unproven choices](#unproven-choices)

## Notes by topic

### Surface operations

Alignment of the control row beside a shared boundary on the line through the boundary pole and the
neighbour's next pole, so the control hulls make one grid; orthogonal projection of pole offsets onto the affine
constraint space; polyspline development of an oriented quad graph in integer UV; face rebuild by interpolation
checked on an independent grid; the number of poles a deformation target needs from the C2 cubic interpolation error
`5 h^4 max|f| / 384` and the end-slope term `h^3 max|f| / 4`; face-separation proofs by Bernstein bounds and, for a
plane beside a surface of revolution about its normal, in closed form (the surface meets the plane in circles; a band
of width `gap / sin t`).

### Alignment input-shape prior

`alignSurface` accepts `inputShapeInfluence` w in [0, 1] (finite, default zero). It acts only when extra blend rows
are requested; zero influence or zero extra rows runs the plain solve with no added matrix entries. The moving
surface and reference are copied and only a certified result is stored.
The construction is Tikhonov regularisation of the free poles about their original positions [Narayan22 slides
6-7]. Let P0 be the pole vector after exact degree and knot refinement, before alignment, and F(P) the bending-plus-
variation quadratic with scalar-coordinate matrix H on the n free poles. The objective is

    F(P) + w a sum_free ||P_i - P0_i||^2,  a = trace(H) / n,

subject to the continuity conditions and fixed remaining poles. Trace normalisation (taken before constraint
rows or the prior are added) and the [0, 1] domain are design choices derived here, not a published
surface-specific weighting law. They make w dimensionless and keep rigid-motion and uniform-scale covariance for
a fixed basis; the pole metric is not invariant under change of basis or arbitrary UV reparameterisation.

For displacement unknowns d from the current iterate Pc, add `w a` to each free coordinate's diagonal and `w a
(P0 - Pc)` to its right-hand side; P0 is held through all nonlinear continuity iterations. Cost is O(N) time and
storage for N total poles, with no new poles and unchanged band; the term is applied after the energy cache so
no prior is cached.

The condition/multiplier solve and final seam check still apply: the prior cannot authorise a seam miss, and a
surface that already meets the continuity is returned as it is. *Derived bounds.* H is positive semidefinite (it
integrates squared spline derivatives), so for `a > 0, w > 0`, `H + w a I` is positive definite on the free
coordinates and a nonempty affine constraint set has a unique minimiser.

For minimisers P1, P2 at `w1 < w2` on the same feasible set, adding their optimality inequalities gives `(w2 -
w1)(D(P2) - D(P1)) <= 0`, D the squared pole displacement: total displacement is nonincreasing, while individual
poles or curvature need not vary monotonically. Nonlinear G2/skew iterations change the linearised feasible set,
so this proves each fixed-row problem, not global convergence or monotonicity there. The final seam
certification is sampled, not a whole-domain Hausdorff proof.

### Collinear alignment: units and rational boundary rows

Collinear alignment uses exact degree elevation and the union of the two normalized boundary knot vectors with
the larger multiplicity at each knot, and edits rows directly (no degree-five fitting restriction). Knot
insertion preserves the parameterized surface [PMC09 §1.4.3, (1.74)-(1.76)];
`Geom_BSplineSurface::IncreaseDegree` preserves the geometry and `InsertUKnot`/`InsertVKnot` with `Add=false`
raise a knot to the requested multiplicity [OCCT8 Geom_BSplineSurface]. *Derived:* for unit row directions a, b
at angle theta, `|a+b| = 2 cos(theta/2)`.

The symmetric fold check compares this dimensionless quantity with 1e-9 (a conditioning choice retained from the
unit-scale algorithm); physical short-row degeneracy is checked separately. Scaling model lengths cannot change
this angular verdict. Rational supports need clamped boundary knots; rows carry weights as well as Euclidean
poles. After exact refinement the moving boundary and next row receive the reference's corresponding weights;
deeper rows keep theirs.

From the rational surface definition [PMC09 §1.5, (1.91)] (derived): write the boundary as `S = A0/W0`; its
cross derivative is `c (A1 - S W1)/W0`, c the positive endpoint degree/knot factor. If the reference next-row
weights are a constant multiple of its boundary weights then `sum N_i w1_i P0_i = S W1`, so placing the moving
next poles at `P0_i + r (P0_i - P1_i)`, common `r > 0`, gives the opposite proportional derivative. This proves
tangency for the constant-ratio option, rational curved boundaries included.

The proportional-weight condition is checked to 64 machine epsilons relative (an arithmetic allowance); a
nonproportional reference weight row with that option is unsupported. Without the option, pole lines align but
neither polynomial nor rational variable row ratios guarantee exact G1; skew changes the target lines and is not
covered by the zero-skew argument. At very large coordinates independently evaluated exact isocurves and surface
points can differ by more than OCCT's default absolute confusion tolerance.

The face construction may then carry the requested positional tolerance onto its new face, edges and vertices if
the default metadata does not validate, and rechecks BRep validity and the interior check. No poles, knots,
weights or curves are fitted or moved by this step and no tolerance larger than the caller's allowance is
introduced; the option is off for other consumers of the face builder.

### Rectangular sheet extension

Extension requires a valid single-wire face whose boundary is its parameter rectangle. `GeomAlgoExt::parameterRectangle`
checks all boundary curves against the four sides: line endpoints, or the complete positive-weight spline pole hull,
which implies the property by the convex hull property applied to each constant coordinate [PMC09 §1.4.2] (derived).
Area agreement cannot establish the condition, since a narrow notch can have arbitrarily small area and finite
depth. Unrecognized boundary representations are not assumed rectangular.

### Unproven choices

Work budgets and halving depths of the separation proofs, the crease threshold, interpolation grid sizes,
and the choice of pole counts ("twenty-one for any target"). The sampled checks are not Hausdorff certificates.
