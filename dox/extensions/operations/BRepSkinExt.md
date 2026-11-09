<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepSkinExt: skins, sweeps and drafts

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Skins and sweeps](#skins-and-sweeps)
- [Unproven choices](#unproven-choices)

## Notes by topic

### Skins and sweeps

Frames along a path are rotation-minimising by the double reflection method [WangJuttler08], with a fourth-order
error in the step; stock `BRepOffsetAPI_MakePipeShell` twists a profile only where the path starts square to it.
A skin's cross-section lines are natural C2 cubics, whose minimum of `integral |s''|^2` over H2 interpolants is
proved in [utilities.md](../utilities.md).

Closed sections are made fit by cyclic seam alignment: the cyclic cross-correlation of two sampled rings is
maximised over start shifts and orientation [Dogan25] (the rotation proof by Kabsch/SVD read; the version here
fixes rotation). Draft of walls that `BRepOffsetAPI_DraftAngle` does not take (it handles planes, cylinders and
cones) uses the ruled isocline surface over the wall's neutral curve, whose every tangent plane makes the draft
angle with the pull direction.

#### Gordon network

Gordon [Gordon71] as presented by [Siggel18 §3.2.1-3.2.3, (8)-(17)]: the profile skin plus the guide skin minus
the intersection-grid tensor surface; also [PieglTiller97 ch. 10, p. 485]. Equation (8) is the compatibility
condition: each skin and the grid must use matching interpolation parameters and bases, which degree elevation
and knot insertion align without changing geometry. Compatibility and rational representability limits of
`GeomFill_Gordon` apply.

A guided skin through three or more profiles first constructs the compatible Gordon network, then applies the
G1/G2 end correction below, avoiding assembly and factoring of a dense plate system for a network already
representable by spline interpolation. The plate construction remains when the network or its correction cannot
be represented or fails its checks. An accepted plate can have a different interior and bending energy, so
validation compares the delivered surfaces; exact curve interpolation alone does not establish equivalence.

Coplanar networks with planar supports keep the plate construction (a plane has zero geometric bending energy,
and rebuilding its unit network rectangle would change an accepted carrier and trims); the model tolerance
defines the coplanarity search and the plate checks decide delivery. For curved networks the accepted interior
is compared with a matched reference surface in both sampled directions.

#### G1/G2 end correction (derived)

Two surfaces sharing an edge and their tangent planes along it meet G2 when their normal curvatures agree in one
direction across the edge (linkage curve theorem [PegnaWolter92]; for the geometric interpretation and its
tangent-plane-continuity assumption see [PMC09 §8.1.3.5]). The skin already holds the profile, so G1 needs only
the way across, `S_v`, turned into the support's tangent plane, and G2 only the normal part of `S_vv` set to the
support's second fundamental form at `S_v`; both moves are along the support's normal.

The correction is `h1(v) m1(u) + h2(v) m2(u)`, h1 and h2 the quintic Hermite functions of the end (`h1' = 1`,
`h2'' = 1` there, the rest zero) that vanish to the second derivative at the nearest place of the next profile.
The other profiles are kept and the skin stays C2; `m1`, `m2` are the moves less their straight interpolation
between the guides, so the guides are kept too. Its sampled field approximation is not a continuous error
certificate.

Where a guide leaves the end profile out of the support's plane, or with another bend, the surface follows the
guide. Nets whose skin carries weights are outside this correction. Exact interpolation does not imply less
bending energy (a counterexample raised the integral from 0.9364 to 1.2673), so this candidate needs a
family-wide geometry and fairness review before it replaces accepted results.

#### Pairing of unlike closed sections (open)

Approximating one conic by non-rational cubic pieces is not the problem of correspondence between unlike
sections. [Floater97 §2, Thm 3.1] give geometric Hermite interpolants and a Hausdorff bound for positive
normalized middle weights up to 3 (for degree three the subdivision order is six); it controls the conversion of
each section, not which points on two sections are joined. Repeated bisection of three initial parameter regions
can explain a knot pattern without determining the correspondence.

Rational reparameterization supplies none either: [HuWangWang18 §2.1, (3)] change weights by `w_i -> lambda^i
w_i` under `t = lambda u / (1 - u + lambda u)`, `lambda > 0`, preserving degree and geometry; their §2.2, (7),
(9), (12) use curvature maxima as breakpoints, normalized arc lengths as new breakpoints and matching
derivatives, as choices for one curve's speed, not a theorem selecting points between unlike curves; a
three-region knot pattern does not establish Möbius pieces (a one-map diagnostic through three places per region
missed the observed conics by 0.234 to 0.550 model units: evidence against that hypothesis, not a bound on
piecewise reparameterizations).

Compatibility is distinct from choosing correspondence: [PieglTiller97 §10.3, pp. 458-459, (10.7)] assume common
section degree and knots and interpolate the control points (homogeneous for rational sections);
[EngleitnerJuttler18 §2-4] permit incompatible section spaces through patchwork splines, with existence in the
chosen space (Thm 3, under shadow-compatibility and smoothness-monotonicity) and a constrained energy
minimization (§4.3). In every case the sections' parameterizations are given, so none identifies an unknown
pairing from an output knot pattern. No pairing law is installed.

#### Solid publication

Planar end caps are the only cap construction of the smooth section skin. A failed cap is a
construction failure for a solid request; an open sheet request returns its sheet. Before publication closure is
checked with `BRepCheck_Shell::Closed`, orientation with `BRepLib::OrientClosedSolid` and validity with
`BRepCheck_Analyzer` [OCCT8 BRepCheck_Shell, BRepLib]; disconnected shells are not silently discarded. These establish
topological closure and validity, not freedom from geometric self-intersection. Nonplanar end caps are open.

### Unproven choices

Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean
fallbacks, the sampling counts of the offset and skin validity checks.
