<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepPatchExt: boundary-first fill

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Problem](#problem)
- [Measures and regularity](#measures-and-regularity)
- [Guides, witnesses and trim checks](#guides-witnesses-and-trim-checks)
- [Energy assembly](#energy-assembly)
- [Open boundaries](#open-boundaries)
- [G2 at a corner](#g2-at-a-corner)
- [Spatial fit and chart](#spatial-fit-and-chart)
- [Stop and ladder](#stop-and-ladder)
- [Reading up to the corners](#reading-up-to-the-corners)
- [Size](#size)
- [Ordering, chart covariance and symmetric requests](#ordering-chart-covariance-and-symmetric-requests)
- [Unproven](#unproven)
- [Stock counterpart](#stock-counterpart)
- [Source files](#source-files)

## Source files

| File | Stage |
| --- | --- |
| `BRepPatchExt_Patch.cxx` | Entry: the sides in turn, the constructions in order, the checks every result passes |
| `BRepPatchExt_Fill.hxx` | What a loop asks, what a fill built is, the stages shared between the files |
| `BRepPatchExt_Loop.cxx` | Mean plane, polygon, and inside tests of the loop seen flat |
| `BRepPatchExt_Net.cxx` | Fit conditions, the least-squares fit, net knots |
| `BRepPatchExt_Measure.cxx` | Misses along edges and at corners, interior and sewn readings |
| `BRepPatchExt_Trim.cxx` | Sides laid into the surface's parameters, trimmed face, spline rim frame |
| `BRepPatchExt_Exact.cxx` | Plane and exact continuations (torus, ball, plane/round/plane, spline face) |
| `BRepPatchExt_Quad.cxx` | Four sides: Coons patch with G1/G2 inner poles; chained quads |
| `BRepPatchExt_OverPlane.cxx` | Any loop: height over the mean plane, trimmed |
| `BRepPatchExt_OverChart.cxx` | Any loop that is a height over no plane: fit over its own chart |
| `BRepPatchExt_Adjacent.cxx` | Adjacent: fit over a conformal chart of the loop and its faces |
| `BRepPatchExt_Chart.cxx`, `BRepPatchExt_Inside.cxx`, `BRepPatchExt_Chain.cxx` | Chart of a loop, inside-fault check, open boundaries |

## Notes by topic

### Problem

Given a closed loop of edges, each optionally asked to meet a neighbouring face to G0, G1 or G2, and optional
interior curves and points, build one face. Candidates are tried in order: a plane; an exact rational
continuation of a torus or ball; a four-sided surface with the sides as its edges (a Coons patch whose inner
poles are moved to meet G1 and G2); for any loop, a fitted surface over the loop's mean plane (or over a chart
of its own, below), trimmed by the loop, with a net refined where the rim turns fast.

The fit minimises a weighted energy of squared second and third parameter derivatives subject to the
constraints. The normal equations are positive semi-definite for a connected net; they are solved by a skyline
Cholesky factorisation with a multiplier outer iteration.

Chart coordinates over the loop are cotangent-weight (Dirichlet-energy) coordinates, injective for a convex
boundary [Floater03, Thm 4.1 (statement read, proof not)]; for a non-convex boundary injectivity is not
guaranteed [Karni, Gotsman and Gortler; Weber and Zorin], so the result is looked at inside its loop before it
is returned. Corner compatibility at a vertex where several patches meet follows the vertex-enclosure condition
[HermannPetersStrotman09].

### Measures and regularity

#### Complete curvature normalization (derived)

A seam compares the two second fundamental forms in one orthonormal
tangent frame. For a symmetric form `[[a,c],[c,b]]` the spectral norm is `|a/2+b/2| + hypot(a/2-b/2,c)`. The
absolute mismatch is the norm of the difference; the relative mismatch is that norm divided by the larger norm of
the two complete forms (capped at one). Only two exactly zero forms have zero normalization: the saddle
`z=k*x*y` at the origin has norm `|k|` although both diagonal entries vanish. This solves the 2x2 symmetric
eigenproblem with the normal-curvature characterization [PMC09 §3.3-3.4, (3.26)-(3.43)]. It is a normalization,
not an approximation-error bound.

#### Parameter units (derived)

Surface-jet rank is tested on independently normalized tangent columns. Their cross
product is the sine of their angle; zero-length columns or a sine at floating-point resolution (eight machine
epsilons, a roundoff policy, not a geometric tolerance) are unresolved. Coordinates use the cross-product dual
basis when the Gram determinant is poorly scaled. The dual identities follow from `x.(y cross (x cross y)) =
|x cross y|^2` for unit columns; the second fundamental form then acts on physical tangent vectors [PMC09 §3.3].
This removes absolute UV-unit cutoffs without asserting regularity of near-dependent columns.

#### Sampled normal rotation (derived diagnostic)

Opposite normals at two separated observations trigger midpoint subdivision, not a fold verdict: a regular
cylinder can turn by more than a right angle between samples under a nonuniform chart (exact rational circle:
[OCCT8 GeomConvert]). Subdivision continues until consecutive normals agree in sign, a sampled differential
loses rank, or UV resolution is exhausted; the last two are reported as unresolved regularity.

This is not a certificate of local regularity everywhere or of global injectivity; the separate Bernstein
regularity certificate on chart fills remains required. Automatic candidate selection prefers a graph when a
four-sided surface has coarse opposed normals; this is a preference, not a fold verdict, and the rotating
candidate is checked with adaptive resolution if the preferred construction cannot deliver.

#### Mean-plane degeneracy (derived)

The centered polygon-area sum is unresolved only below `8*n*epsilon*perimeter^2`
(n the number of vertices), not below a fixed share of squared perimeter, which would impose an aspect-ratio limit
on exact planar rectangles. Each centered radius is bounded by the perimeter, so bounding the cross products and
accumulated additions gives this conservative first-order roundoff policy. It is no lower bound on geometric area
or on the existence of a fill; the later geometric and tolerance checks remain required.

#### Planar collars and graph directions (open construction; derived)

A planar support with unit normal `n` is a height over a plane with unit normal `d` exactly when `n.d != 0`.
Orthogonal projection restricted to its tangent plane has singular values `1` and `|n.d|`, so its inverse
condition number is `1/|n.d|` and the height gradient has norm `sqrt(1/(n.d)^2 - 1)`. The `0.8` graph-selection
threshold is therefore a conditioning policy, not an existence condition. Three mutually perpendicular support
normals have best common margin `1/sqrt(3)`.

A simple projected boundary and a global single-valued interior must be established separately. For a polynomial
tensor graph, affine reproduction gives exact planar boundary jets when every height coefficient whose closed
basis support meets a planar rim is set to that plane's affine height at the tensor Greville abscissae: local
support makes the value and first two derivatives planar at the rim, degree five with simple knots gives C2.

Incompatible planes must not prescribe the same coefficient; refinement separates disjoint collars but cannot
resolve inconsistent jets at a shared corner. The identities follow from [PMC09 §1.4.1, (1.58)-(1.61)] by
polynomial reproduction in both tensor directions: with `xi_i = (t_(i+1)+...+t_(i+p))/p`, `p
(xi_i-xi_(i-1))/(t_(i+p)-t_i) = 1`, so by (1.61) `(sum xi_i N_i,p)' = sum N_i,p-1 = 1`; the clamped left
endpoint fixes the constant, `sum xi_i N_i,p = u`, and partition of unity reproduces any affine height.

The construction is not enabled: in a three-direction planar-rim experiment (4 761 to 68 121 coefficients met
the rim checks) the interior-curvature ratios were about 15, 11 and 10 against the allowance 4. These are failed
candidates, not a proof of infeasibility; exact planar jets do not establish interior fairness, topology or
agreement with the intended shape.

### Guides, witnesses and trim checks

#### Guides with positional tolerance (derived)

A guide is a positional constraint, not an exact endpoint two-jet, so an endpoint tangent or normal-curvature
mismatch cannot veto construction: `(x,y,z)=(.01*t,50,h*t)`, `0<=t<=1`, lies within `h` of the planar square for
every `h<=tolerance` whatever its endpoint angle. Construction and the final measurements decide feasibility at
the requested tolerance. This asserts no existence for arbitrary guides. Endpoint jet information is kept only
as a hint explaining a fit that already failed its result checks; it never rejects a delivered fit.

A rejected adjacent fit marks angular and complete-form curvature misses as well as positional ones, and the
between-row side checks call the same measurement and easing alternatives as acceptance (derived). Finite
refinement does not prove existence or convergence for arbitrary boundary jets.

#### Owned witnesses near curved trims (derived)

A fixed-resolution trim polygon can misclassify points near a curved rim. An interior-distance witness uses the
parity of boundary arcs certified to avoid the query point: positive-weight Bezier hulls enclose each arc in a
convex tube about its chord [PMC09 §1.4.2, convex hull property]; if the tube excludes the point, straightening
the arc to its chord cannot cross the point and preserves crossing parity [HormannAgathos01 §1-2, (1)-(3)].

Only uncertain arcs are subdivided; uncertainty is kept once hull deviation reaches the floating-point enclosure
guard. A parameter rectangle is a witness only if its centre is inside the delivered trim and no boundary
enclosure intersects it. This is not interval arithmetic; unresolved conversions, offset approximations and
near-boundary points give no witness. Adjacent owned pcurves can represent their common vertex differently.

A parity chain adds connecting segments only where the gap fits the enclosure guard; for larger gaps ownership
comes from `BRepClass_FaceClassifier` on the delivered face, not an invented closing curve. Before ownership is
propagated over a nonzero rectangle, every boundary-arc tube and the whole coordinate box of every endpoint gap,
enlarged only by the guard, is excluded; unsupported pcurves cannot supply a cell certificate, and points
classified on the rim stay unresolved.

The propagation is derived: an open connected rectangle disjoint from the boundary lies in one component of its
complement, so one authoritative interior point suffices. It assumes the face's own ownership semantics and
certifies no malformed wire. It changes neither the cells selected for surface regularity checks nor proves
injectivity of a fitted face.

#### Interior trim checks (derived)

Requested points are measured against the delivered face, boundary included. Positive-weight guide splines are
subdivided until their control hull has a tolerance-ball cover about a measured face witness: with `d` the
witness distance and `r` the hull radius, `d+r` bounds every point by the triangle inequality; for a convex
planar trim the tolerance neighbourhood is convex, so every hull pole suffices [PMC09 §1.4.2].

For polynomial guide and surface pieces a paired cover avoids tolerance-sized balls: surface cells are
restricted to an in-trim UV rectangle and its diagonal, with the certificate partitioned at the guide's knots
and at every crossing of a surface knot by the UV chord; exact knot insertion and restriction give polynomial
intervals [PMC09 §1.4.3, (1.74)-(1.76)]. For bidegree `(p,q)` the diagonal poles are `R_k = sum_(i+j=k) C(p,i)
C(q,j) P_ij / C(p+q,k)` [PMC09 §1.3.6, (1.56)].

After degree elevation the largest distance between paired guide and trace poles bounds the whole difference
curve by the convex-hull property. The bounds are exact-arithmetic arguments; floating-point conversions and
pole calculations are not interval certificates. A local projection only supplies an in-trim upper-distance
witness; failing that, the global trimmed-face distance query is used. For non-spline curves requiring
approximation, the converter's reported maximum error is reserved in the budget (not an independently validated
bound). Unresolved coverage is refused, never accepted from finite guide samples; coverage can be costly on
curved faces at tight tolerances.

#### Reusing projections (derived)

A projection onto a construction's carrier is not an owned witness on the delivered
trim: it can seed a UV search, but the surface value and trim ownership must be checked. Finite samples do not
bound a guide between samples; positive clearance from the trim establishes ownership of an enclosed UV path, not
the spatial distance from the guide to it. The final cover may be skipped only given an already certified
whole-guide cover on the same delivered surface and trim, or a proof that every later operation preserves it.

#### Local extraction

Knot insertion preserves the polynomial image in exact arithmetic [PMC09 §1.4.3]. Both
constructors of [OCCT8 GeomConvert_BSplineSurfaceToBezierSurface] copy the input surface before restriction and
knot insertion, so a bounded conversion per visited cell is not guaranteed cheaper than one whole conversion.
Ownership tests and local projections, not conversion, dominated the cover cost of a guided polynomial panel;
neither projection-seed reuse nor lazy bounded extraction established cost neutrality (measurements, not
complexity bounds).

### Energy assembly

#### Assembly (derived)

A quadrature cell touches local poles `(a, b)` in order `p = (firstU + a) countV + firstV + b`; these indices
increase strictly, so an element's lower triangle lies in the global lower triangle. Each given pole's
coordinate (or its dot product with the height direction) is read once per fit and each element's
local-to-global indices once.

The energy does not couple spatial coordinates, so visiting one coordinate at a time keeps each right-hand-side
entry's sequence of subtractions, and each matrix entry receives the same additions in the same order (zero
entries stay omitted from the skyline). Only indexing and repeated projections are removed; quadrature,
reductions, constraints and refinement decisions are unchanged. Four fixed corners determine the affine
nullspace, so the zero-energy affine plane is recovered from perturbed poles.

### Open boundaries

An open boundary asked for Quad sided, whose runs and gaps make four sides with corners, is filled by the
four-sided surface untrimmed: the poles of a gap side are free like the inner ones and nothing is asked along it
(a side with more poles at degree 5 than a fitted net may have is first brought down to degree 5 within a tenth
of the tolerance, as an arc is). If that quad is not a face to keep, the boundary is filled as any open one is.

Where the two rims face each other across two gaps, the fill must turn from the tangent plane of one to that of
the other, often more tightly than either rim bends (the end of a tube onto the end of a plate). The inside is
then also judged against the *bridge* between the rims: along each crossing, the Hermite curve that leaves and
meets each rim square to it in its tangent plane, at the speed of the chord, with the face's bend that way where
G2 is asked.

It is built from the rims' places, tangent planes and bends only, so it is pose- and scale-covariant, and a
bulge, ripple or spike that leaves the rims unchanged leaves it unchanged. The fill may bend four times as
tightly as the rim or the bridge within a tenth of the crossing's length of the same place (both factors
judgements). With no gap, or non-opposite gaps, the reference is the rim alone.

#### A section as reference

After a guided height fit fails, a planar section guide between two disconnected rims supplies a reference
candidate. Its plane meets each rim at an anchor; quintic Hermite segments extend the guide to the anchors. At a
smooth rim the extension's first and second derivatives come from the support jet, differentiated subject to the
section-plane constraint; an unsupported rim asks for position alone.

The reference is the single tensor spline `S(u,v) = R(u) + C(v) - R(u*)`, adding the two curves' poles and
subtracting the anchor, without a pole-count refusal. *Derivation and scope.* This translation member of
curve-network blending contains the section algebraically: if the support is `R(u) + Q(v)` with its transverse
two-jet independent of `u`, matching the section's two-jet to `Q` at the anchor gives the same surface two-jet
along the whole supported rim.

A rail whose projection derivative on the section-plane normal stays away from zero and a simple regular planar
section give a regular injective translation surface; these are sufficient conditions for this class, not an
existence theorem for arbitrary networks. Opposite rims with independent parameterisations are measured against
their original curves. The candidate must pass the delivered-rim, support, guide, validity, interior and sewing
checks; a failed reference keeps the refusal.

The curve-network Boolean sum and its conversion to one tensor spline are in [Siggel18 §3.2.1-3.2.2, (10),
(14)]; support-derived second-order ribbons and their continuity conditions in [SalviVarady14 §3.3.2, Appendix].
The section-plane acceleration correction and this restricted route are derived here. The interior curvature
allowance is the unproven judgement below; a compatible guide can still be refused.

### G2 at a corner

#### Corner compatibility (derived)

Take a common oriented tangent plane at a corner p. For a regular C2 fill F,
the second fundamental form II_F is continuous there. G2 contact to support A along one incident edge gives
II_F = II_A on that edge, hence II_F(p) = II_A(p); the other edge likewise requires II_F(p) = II_B(p). So
II_A(p) = II_B(p) is necessary: a nonzero difference rules out exact G2 at that corner, however many poles or
internally G2 patches the fill has. This is the elementary corner obstruction, not a sufficient existence theorem.
II and its relation to normal curvature: [PMC09 §3.3, (3.26)-(3.30), Thm 3.3.1]. The scalar twist residual alone
cannot detect it: in a common orthonormal frame II_A = diag(k, 0) and II_B = 0 have zero mixed term but differ by
|k| in normal curvature. This covers a bending strip meeting a flat strip and rounds of unequal curvature forms at
a corner. These are impossible exact requests, not missing capabilities. A seam of lower continuity or a singular
vertex changes the regular G2 request. The proof does not make the relaxed (eased) problem below impossible: failure
to meet the limits beyond the ease is a finite-approximation issue.

#### Ease at a corner

Where two sides asked for G2 meet at a corner whose faces bend differently (an arm of a round meeting a strip
whose edge bends across the arm's axis), no surface twice differentiable at the corner has both bends there.

It has the one bend that the two edges and the tangent planes along them leave it: four numbers (each edge's
normal curvature and the turn of its tangent plane along it) for the three of one bend, solved by least squares,
with each face's bend across its edge weighted a thousandth as much, so that it decides only what the four leave
open (for an edge that runs straight through the corner, the bend across becomes the mean of the faces').

The bend each side asks across itself eases from that one at the corner to its own face's as `1 - (10 t^3 - 15
t^4 + 6 t^5)` of the difference over the ease length: the quintic whose first two derivatives vanish at both
ends, so the bend asked stays twice differentiable along the side and agrees with the other side's to second
order at the corner. Within the ease a miss is allowed a share of no less than the ease itself (near a corner of
straight edges both the bend asked and the surface's are near flat).

The fill is held to the eased bends; the per-edge curvature it reports is against the faces themselves, so the
corner's miss is reported and the fill is not refused. *Ease length (derived).* If the faces' bends are met on
both sides beyond the ease, the fill's bend has to turn across the corner by the whole difference at every
distance d from it (on one edge it is that edge's curvature, on the other the other face's), so it changes by
that over an arc of about d; a net follows it only where d is no shorter than its span.

The ease is one span of the coarsest net (the loop's size over 6), at most half its side (a judgement). Where no
fill that eases so meets the faces beyond it (three rounds round a box corner at tolerance 1e-6; unequal rounds
and their planes round a spherical cut), the refusal names the corner and says G1 can be had. The G2 analogue of
the vertex-enclosure condition is stated open in [HermannPetersStrotman09]; this is the elementary part of it
that a fill can be held to.

G2 is then met on the middle of each side only; coarse nets that meet an eased loop can be lopsided; and the
ease of one coarse span is a judgement bounded by measurement, not an optimum. A flat strip meeting a bending
one is held beyond its ease to the flat G2 limit, which no eased fill meets, and is refused naming the corner.

### Spatial fit and chart

#### Over a chart of its own

A fill over the mean plane is a height over it, with slope `tan(acos |n.z|)` where the rim's faces have normal
n; the plane route refuses a loop where some `|n.z| < 0.17` or where the loop seen along z crosses itself. Some
loops are a height over no plane: where the face normals along the rim sweep a closed half circle or more (a
wing nose, a window round more than half a cylinder), every direction d is square to one of them by the
intermediate value theorem, so every height over every plane is vertical somewhere on its rim.

The nearest plane is read exactly: with the signs s of n.d fixed near the best direction of a search over the
half sphere (a Fibonacci lattice of 4096 directions in the loop's frame, so it turns with the loop), `max_d min
s n.d` is the distance from the origin to the convex hull of the signed normals [Boyd04 §8.1.1, support
function], found by Gilbert's iteration [Gilbert66]; every answer is a direction's own least, hence a bound from
below.

A loop refused by the plane route whose margin, read at 65 places of every side asked G1 or G2 with their ends,
is below the same 0.17 is filled over a chart (`overChart` in `BRepPatchExt_OverChart.cxx`, chart in
`BRepPatchExt_Chart.cxx`); any other loop is filled over the plane, and a loop the chart does not fill is
refused as over the plane.

The loop is first spanned by a disc of triangles: an even lattice in the unit disc, the loop on its rim at even
steps of its length (160 places), whose places in space minimise the sum of squared umbrella Laplacians at every
node of the rim and inside it, the rim held (discrete thin-plate fairing [Kobbelt97; BotschKobbelt04]).

Along a side to be tangent to a face, a ring of nodes outside the disc stands in the face at each piece's
middle, less the way the fill leaves the side (away from the face, read in the face's own parameters) times the
height of an even triangle; with the rim's umbrella near nothing, the nodes inside next to it stand where the
face's tangent plane goes on, the discrete form of the first-derivative condition that two fixed rings give for
the biharmonic equation [BotschKobbelt04]. Elsewhere the outer node is free (a natural end).

The umbrella is the Laplacian of the lattice's own parameters, so the disc is laid out again, twice, as an even
lattice in the chart it lays flat to, where an even triangle is near even on the surface and the umbrella near
its cotangent Laplacian [Pinkall93], and faired again; from the first of these on it is pinned to 33 places of
every guide where it was nearest them (each pin counting 100 umbrellas, a judgement).

It is laid flat by Boundary First Flattening with the scale held at none along the rim [SawhneyCrane17 §5-6.1;
SpringbornSchroederPinkall08 App. E]: angles kept everywhere and lengths along the rim, so that the thin-plate
energy of the parameters is near the surface's bending energy, which it is not over a plane where the surface
stands steep. The flattening has no injectivity guarantee for a non-convex boundary [Floater03 Thm 4.1;
SawhneyCrane17 §7]; a triangle laid flat folded over refuses the chart.

A net free in space is then fitted over the chart as over the plane, each side at its places in the chart, its
tangent plane asked as `N.S_u = N.S_v = 0` and its bend as `n.S_nn = II(S_n, S_n)` about the surface of the pass
before (both conditions on derivatives in space, none on how the chart lies), and trimmed by the sides' curves
in its parameters. Such a net can meet every tangent plane by not leaving the rim at all (`S_n` running to
nothing), which is no parametrisation.

In a chart that keeps the rim's lengths, a surface to which the chart is conformal leaves the rim as fast as it
runs along it, `|S_n| = |S_t| = 1`; along each side asked G1 or G2 the speed along the way the fill leaves is
held within 0.5 of that, asked only where it leaves that band (an inequality, not a wish; the band is a
judgement, since the chart is the disc's, not the fill's).

The net starts as the coarsest of the ladder and spans are halved where a place is missed, as for the Adjacent
fit, up to the poles a net free in space is let have (`freePoles`, 640); once a net meets the limits, the next
finer is fitted and the coarser kept where they lie within `settled` of each other (the stop below).

The surface is kept only where `S_u x S_v` is shown away from zero on every cell used: per cell it is a
polynomial of degrees (2p-1, 2q-1) whose Bernstein coefficients come from those of `S_u` and `S_v` [Farin02
§5.8]; by the convex hull property a cell whose coefficients all have a positive part along one direction has no
zero, and a cell that does not is cut in four, at most four times (the surface form of the Bernstein
Jacobian-sign test [Gravesen14]).

#### Cost forecast (derived)

`netSpacing` predicts a chart's cost by replacing higher derivatives by those of the osculating circle: a circle
of curvature k has arc-length derivative magnitude `k^p` at order p+1. This is exact for the circle.

Curvature at a point of a general curve or surface bounds neither its higher derivatives nor those of its chart,
so the forecast and its 640-pole cutoff are computation policy: no approximation-error theorem, no necessary
pole count, no proof that a refused fill does not exist. [SandeManniSpeleers20 §1.2, (2)-(3), §6-7] give upper
error estimates in the actual Sobolev seminorm and, for mapped geometry, derivatives of the map; those are
L2/Ritz estimates, not pointwise curvature bounds, and cannot be inverted into a lower bound on a feasible pole
count. Replacing the surrogate requires measured constructions and the same final geometry checks.

#### Spatial curvature condition (open)

An unconstrained spatial fit enforces one transverse curvature equation, with five outer linearizations. G0 and
G1 tolerance do not bound the omitted tangent and mixed second fundamental-form components. A complete
correction must check the actual nonlinear form and the delivered parameter curves, including places between
fitting samples, and establish regularity and interior acceptance. Support-derived second-order ribbons give
local jets [SalviVarady14 §3.3.2], not a convergence theorem for one fitted spline.

Globalized constrained iterations and feasibility restoration are in [NocedalWright06 §18.4, Alg. 18.3]; their
hypotheses and convergence must be verified for the assembled problem, and merely adding three frozen curvature
rows establishes neither. Open: a constrained feasibility solve with complete curvature conditions and measured
convergence, then refinement where the verified geometric residual requires it.

The cap of a cone standing 80 degrees or more against its circle asks for 529 to 576 poles at tolerance 1e-4,
within the 640; its coarse nets are not shown regular, and it is refused. The chart's own map is regular (a
folded triangle refuses it); what is not shown regular is `S_u x S_v` of the net fitted over it, which the chart
does not decide.

#### Chart constants not derived (judgements)

160 places round the loop, 33 pins per guide, two relays, the pin's
100, the band's 0.5, four cuts of a cell, the 4096 directions (they only seed the exact step), and the 256 bins over
which the spans asked are read across the chart; the Favard constant is taken at its upper end, so that the spans
are, as far as the bound goes, enough.

### Stop and ladder

A net that meets every limit is not yet the fill asked for: that is the fairest surface of all that meet them,
which the nets of the ladder approach as they become finer.

Once a net meets the limits the ladder continues until two nets in a row lie within `settled = 1e-3` of the
loop's size of each other, read from twice as many places as the finer net has spans each way inside the loop to
the nearest place of the other surface; the coarser of the two that meets the limits is kept, and where the
ladder or the allowed nets end first, the finest that meets them.

This is a Cauchy test: where the nets converge at least twice as fast as the ladder refines (ratio `rho <=
1/2`), the kept net is within `d / (1 - rho) <= 2 d` of their limit, d the distance between the two. A share of
the size makes the stop the same at every size and pose.

A net free in space and a height over the plane minimise different energies and approach different surfaces;
they are not one sequence, and once nets free in space meet the limits no height is tried after them. `settled`
is a judgement (about the height by which a tilt of 0.1 degree moves a surface over half the loop's size). The
cost is at least one net beyond the one kept; nets that a coarse net met only roughly come out finer. Where nets
keep differing by more than `settled` the whole ladder is run.

Open: the Cauchy bound needs the distances to fall, and they need not: on some loops (an eight-edge G0 opening
with nets 2.38e-3, 7.0e-4 and 7.5e-4 of the size apart) the ratio is above 1/2, the nets stand on a plateau
below `settled`, and the stop is a judgement, not a bound.

#### A ladder that does not settle

The distances between nets may not shrink at all (nets 3.8e-3, 3.1e-3, 1.3e-3, 2.9e-3, 1.2e-3 of the size apart,
up to the finest). A condition is met once it is within its reach: the multipliers stop at the first shifts that
put every condition there, and the correction then moves the net onto the curves by the least movement of its
poles. What a net delivers is one surface of the band of all that meet the limits, at which that iteration and
correction stopped, not the least bending of the band.

The band has an inside (the limits are inequalities) and as nets become finer nothing draws their surfaces to
one point of it. What is asked is the fairest surface that meets the limits, the least of the bending energy
`integral (k1^2 + k2^2) dA` over the band. So where no two nets in a row come within `settled`, of the nets that
met every limit the one of least bending is kept (`bendingInside`, read at the places the inside of a face is
looked at; `k^2 dA` has no unit, so the choice is the same at every size and pose to rounding).

A ladder that settles keeps what it kept. The energy is the measure of fairness and the places are those of the
inside check; nothing is newly tuned. Pose dependence is measured, not removed: no change in whether a case is
refused was seen over rigid-motion variants.

### Reading up to the corners

A fill is measured at the 129 places of Measure Continuity along each side and, between the first and last of them
and the side's ends, at places that halve their way to each end down to `cornerInset` (1e-6 of the side; a fraction,
so the same at every size; a judgement of how near is at the corner). G0 and G1 are read there as everywhere; the
bend is read there too, but not within an ease the corner's faces are shown to ask for (see `agreeAtCorners`),
where it is read from the first place on. Fills that pass the 129 places can leave their rim near the corner.

Where two sides that are only to be met meet at a corner p, a surface holding both is tangent at p to the plane
of their two tangents: each edge lies on it, so each tangent lies in its tangent plane, and two non-parallel
tangents span it. The plane may stand far from the surface around it, but it is asked at p only: along each edge
the surface's tangent plane need hold only that edge's tangent, and turns as fast as it will. Where the edges
run on within the G1 angle the two tangents span no plane and nothing is asked.

Within the tolerance a surface need not be tangent to it, and asking it of every such corner bends fills that
were within before; so a correction onto the sides holds, besides its even places, those of the measure's places
near such a corner that the corrected surface misses by more than the tolerance, and finer nets of the ladder
are asked the plane at that corner.

### Size

Lengths are shares of the loop's size or of the tolerance, with one exception by design: where a fill goes over
into planes and other faces at once (or into a cone alone) its rim and its edges' parameter curves are made to a
tenth of `Precision::Confusion()`, OpenCascade's absolute least distance, so that its edges keep that as their
tolerance and STEP gives the joint back. A fitted net moved onto its curves ends its ladder when the bounded rim
correction leaves its own places further off than the tolerance (the correction's own residual, no constant).

The orientation of two collinear pieces of one straight side is decided by the sign of the `orient2d`
determinant only where it exceeds the bound on what rounding its products and their difference can do
[Shewchuk97, ccwerrboundA], a bound relative to the numbers themselves and so the same at every size.

At a corner where both edges run straight on rounds eased to one bend, everything the corner's four numbers ask
is rounding; their share of disagreement counts only above the rounding of the faces' own curvatures, `2^-40` of
the larger principal curvature (`jetRounding`, a judgement: a few thousand units in the last place, far below
any curvature asked).

The margin is thin: rounding alone came to at most 9e-17 of the larger curvature against the floor 9.1e-13, but
the least genuine disagreement met, 6.7e-12, is only 7.4 times above it, so a genuine disagreement eight times
smaller would be taken for rounding and its turns dropped.

A floor derived from the faces' own representation (each derivative off by up to `(gamma_n + eps) sum |N_i^(a)| |N_j^(b)|
|P_ij|` over the poles where the face stands, carried into the normal curvature) was measured six to
ten times above `2^-40` on flat strips and is not used: strips interpolated from posed points are rounded by the
solve's condition, not by the poles' storage, and no bound from their representation tells that.

A floor at the flat share of the G2 limit would be the measure's own and is the structural path once the fold of
the net over the plane without those turns is fixed. Open: a saved opening of near-flat spline faces (curvature
6e-7 on a loop of 38) is filled only when the turns asked by its small but genuine disagreement are kept;
without them its fill folds over, an open defect of the net over the plane, not of the corner.

At 1e-3 size the corner of three rounds (tolerance below `Precision::Confusion()`) is made and one
unequal-rounds corner is refused at 1e3. The switch from the height graph to the spatial net at a support normal
of 0.8 against the mean plane changes a cone cap's area by 2.1 % between 36.8 and 36.95 degrees, and an analytic
cone's cap changes between two shapes from 70 to 74 degrees; a cap over a circle on a cone is refused as too
steep from about 80 degrees.

### Ordering, chart covariance and symmetric requests

Permutation invariance of an annotated boundary differs from isometry covariance of its fill, `F(g X) = g F(X)`,
with supports, guides and continuity orders transformed together.

If `g X = X` a covariant construction preserves that symmetry; a unique minimizer supplies this only when its
feasible space and complete objective are also preserved (convergence of a penalized local solve does not prove
uniqueness). *Frames (derived).* A unique directed frame cannot be covariantly selected at every symmetric
request: if `h X = X` and `Q(h X) = h Q(X)`, an invertible frame gives `Q(X) = h Q(X)`, forcing `h` to be the
identity.

A chart defined up to an orthogonal coordinate change works when all later stages respect that change.
[DymLawrenceSiegel24 §2.2-3] discuss further continuity obstructions and weighted frames. The proxy chart is not
such a complete construction: changing its seam rotates its circular boundary against a regenerated triangular
lattice, and a mesh whose connectivity respects cyclic shifts and reversal still cannot resolve an ambiguous
guide preimage (one of two equally close triangles can break a symmetry).

A concentric-ring mesh with four triangles per annular quadrilateral removes the arbitrary diagonal but is valid
only while its homothetic rings do not fold. An experimental guide initializer averages triangle-foot chart
coordinates with weights `area * exp(-(distanceSquared - minimumSquared)/meanArea)`; normalized weights are
invariant under isometries and uniform scale and the average transforms under affine chart changes (derived).

The mean-area bandwidth is an empirical choice, not an error bound; the average need not be a nearest preimage
nor lie inside a non-convex chart. Neither is a production caller. [SawhneyCrane17 §4.2, 5] give the flattening,
not a guarantee for these mesh or correspondence choices. A fixed tensor spline space is generally not closed
under chart reflections: reflecting `u^5 v^5` about an axis at pi/8 gives `(u^2-v^2)^5/32`, outside bidegree
(5,5). Coordinate reversal requires reversal-symmetric knot vectors, which arbitrary adaptive knots do not
ensure.

Complete total-degree disk polynomial spaces, including complete real Zernike levels, are closed under
orthogonal chart changes [GreengardSerkh18 §2.3, 4]. A separate spline approximation needs its own budget:
approximating a reflection-symmetric surface within Hausdorff error e implies mirror distance at most 2e by the
triangle inequality (derived), with no derivative or regularity bound.

### Unproven

The interior acceptance gate (a fault when the interior bends more than 4 times the rim's bend, with a floor of
2/size) has no derivation; a bound of the form `K sqrt(2D)` in the energy gap D was tried and rejected. The bridge only
replaces the reference that the 4 multiplies; the 4 and the tenth are judgements. Against the bridge, defects that
leave the rims unchanged are caught less often where an elbow turns than against the rims alone (gentler bulges and
spikes pass more often). The flat G2 limit, the rules by which a ladder that has met nothing yet gives up
(`stillFalling`, `same`, `unruly`, which compare misses of unconverged nets and so depend on rounding and
therefore on pose), the span and budget constants of the net, the `settled` share and the 0.17 slope limit are
unproven.

### Stock counterpart

`BRepFill_Filling` and `GeomPlate_BuildPlateSurface`: a plate fit that accepts G0, G1 and G2 point and
curve constraints and minimises a different energy.
