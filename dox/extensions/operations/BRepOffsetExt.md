<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepOffsetExt: offsets, shells and thickening

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Offsets and shells](#offsets-and-shells)
- [Face offset that keeps the topology](#face-offset-that-keeps-the-topology)
- [Local focal witnesses for a sheet's normal sweep](#local-focal-witnesses-for-a-sheets-normal-sweep)
- [Elementary focal distances](#elementary-focal-distances)
- [Regularity of network surfaces](#regularity-of-network-surfaces)
- [A circular tube's parallel body after its hole closes](#a-circular-tubes-parallel-body-after-its-hole-closes)
- [Offsets of moved faces: domains, seams, pipes](#offsets-of-moved-faces-domains-seams-pipes)
- [Clearance coverage and rims](#clearance-coverage-and-rims)
- [Neighbour-bounded planar push-pull](#neighbour-bounded-planar-push-pull)
- [Face-offset endpoint certificate limits](#face-offset-endpoint-certificate-limits)
- [Separated tangent junctions](#separated-tangent-junctions)
- [Transported contact uncertainty is confined to its measured domain](#transported-contact-uncertainty-is-confined-to-its-measured-domain)

## Notes by topic

### Offsets and shells

A solid offset first classifies the faces the distance would shrink to nothing, removes them from a working copy by
extending their neighbours (`BRepAlgoAPI_Defeaturing`) and offsets that copy; a face too tightly bent in part is
refused by name. Stock `BRepOffset_MakeOffset` removes vanishing faces only where every face is planar
(`myIsPlanar`); a fillet whose radius the distance reaches comes back as an invalid solid or not at all.
Push-pull follows [ZouFeng19]: the solid swept by the moved faces is built and joined to or cut from the body by a
Boolean, which resolves whatever the move does to the rest. Thicken divides a C0 sheet's supports without changing
their geometry before offsetting.

### Face offset that keeps the topology

`BRepOffsetExt_Kept.cxx`, tried first by push-pull where every face moves. It holds for any set of faces with no
tangent neighbour that stays; a move of only some faces keeps the engine's route. The face offset keeps every face
(the intersection join). Each face goes onto its offset surface, which `BRepOffset::Surface` gives in the surface's
own parameters: a plane translated; a cylinder, cone, sphere or torus with the same axes and `(u, v)`, since
`S(u, v) + d n(u, v)` is again such a surface in those parameters; otherwise `Geom_OffsetSurface`.

#### Tangent edges

If two faces have the same normal n along an edge c and move by the same d, the point `c + d n` lies on both
offset surfaces with normal n on both, so they touch along the offset of the edge and its curves on the faces
stay. This is tested, not assumed: the two faces' points of each edge parameter after the move agree within the
tolerance, or within the edge's own tolerance where that is larger. Two faces tangent only to an angle theta
come apart by `|d| 2 sin(theta/2) <= |d| theta` more when both move by d.

An input G1 only to theta (its edge tolerance that wide) has its edge kept, and the result's edge carries the
input's tolerance plus `|d| theta_max`. `theta_max` is found, not bounded (judgement): 8 samples in each piece
of the edge where its curves on the faces are polynomial and cross no knot line of either surface, each local
maximum refined by golden-section search.

A bound would need the faces' curvatures bounded along the edge, and the hull bound on curvature was three
orders loose on blend sheets; sampling found the larger maximum than a 4 097-sample reference on every nearly
tangent edge tested, while 33 even samples fell short by up to 36 per cent. The intended budget is that edge's
carried tolerance and the input vertex's plus the same spread (the endpoint gap limit below is not covered); the
result is never claimed more accurate than its input.

#### Sharp edges

The new edge is the crossing of the two offset surfaces, solved in the plane through c(t) normal to the old
edge: four equations `S1(u1, v1) = S2(u2, v2)`, `(S1 - c(t)).c'(t) = 0`, by Newton from the old edge's curves on
the faces, so it keeps the old parameter t.

The curve and the curves on the faces are cubic interpolants of those sections, refined locally (a midpoint is
added where an interpolant is further than a tenth of the tolerance from the section there); refinement stops at
16 halvings of the first spacing or 4096 samples (judgement: a curve a cubic cannot follow there has a kink).

An offset surface is one order less smooth than its basis (its normal takes a derivative), so where the edge
crosses a line along which a moved surface is less than C2, or the old edge (whose tangent sets the section
planes) is less than C3, the curve is interpolated apart on each side and joined C0 there; the crossing is found
among samples and by bisection.

Newton steps that do not reduce the residual are halved; where two surfaces only touch and the residual stalls,
a point within half the tolerance of both is accepted (judgement: the other half is left to the interpolation, a
tenth, and the input's own edges).

#### Vertices

On a tangent edge's (or seam's) moved curve, crossed with the face around the vertex that meets it most
steeply among those that meet another at a sharp edge there (all tangent: the vertex moves along the common normal);
between two edges of the same two faces, in the section of one of them; else where three faces cross (the three whose
normals are furthest from a common plane). A vertex then not on every face around it is declined by this local
construction; that does not prove the topology must change, since approximate contact continuations can also leave the
other support.

#### Certificate

apart from the construction: a valid solid; points inside every face and along every edge at the moved
distance from their sources' unbounded surfaces, on the side moved to (foot point by Newton, Gauss-Newton where the
Hessian is not positive definite); faces clear of each other (below); volume grown by an outward move and shrunk by an
inward one; and, where every face moves, nothing nearer the body than the distance. A failure of any falls back to the
engine's routes. At a pole, where `S_u x S_v` vanishes against the derivatives' size, the certificate checks the
distance alone (judgement: below a millionth of `|S_u|^2 + |S_v|^2`), as the normal's sign is not defined there.

#### Faces clear of each other

The general self-interference check (`BOPAlgo_ArgumentAnalyzer`) on the exact faces, on whatever the separation
pre-pass (`BRepSurfaceExt_FaceSeparation.hxx`, pole hulls and Bernstein bounds) cannot prove clear. The pre-pass
does not read offset surfaces, so it runs on proxies: a copy whose offset faces are splines within `eps` = a
hundredth of the tolerance of them, each face's tolerance grown by `eps`.

Two proxies proved further apart than their tolerances have their exact faces apart too; that a proxy is one to
one says nothing of its exact face, so each offset face also goes to the check by itself unless the pre-pass
leaves it in a part that is checked. The result keeps its exact faces.

Skipping the check by the reach [Federer59 Thm 4.8] (`|d| kappa_max < 1` per face, faces sharing no vertex
further apart than both moved faces reach) was dropped: a rigorous kappa bound per rational Bezier patch was
sound but three orders above the sampled curvature on blend sheets, and the reach leaves faces that share a
vertex unproven.

#### Offsets that load again

The binary BREP reader builds an offset surface with OpenCascade's continuity check (`BinTools_SurfaceSet`),
which refuses a basis that is C0 at a knot unless it finds it G1 there, so a body holding the offset of a blend
sheet whose parts are joined C0 in their parameters cannot be loaded from a saved file.

Such a face's offset is built unchecked for the construction and, once the face's extent is known, replaced by a
spline over its box within a tenth of the tolerance (`GeomConvert_ApproxSurface`), its edges' curves, vertices'
points and continuity records moved onto it and their tolerances grown by the error.

An edge whose old curve is a line or a circle keeps one: the offset of a line on a plane or cylinder is a line,
of a circle round a surface of revolution's axis a coaxial circle, with the parameter kept. *Read back as
written.* A body read from a saved file has the directions it was written with, bit for bit (`gp_Dir::Restored`
in the BinTools and GeomTools readers; [corrections.md](../corrections.md#gp_dir-restored)).

Open: the X and Y directions of an axis placement (`gp_Ax2`, `gp_Ax3`, `gp_Ax22d`) are rebuilt from the main and
X directions by cross products when read and come back exact only where those products are; a reader that keeps
orthonormal axes as written is not made.

### Local focal witnesses for a sheet's normal sweep

`BRepOffsetExt_Focal` encloses non-periodic rational B-spline and Bezier supports directly from stored poles,
positive weights and knots, on each nonempty knot rectangle. The homogeneous first and second derivative control
nets use the derivative formula [PMC09 §1.4.2, (1.65)-(1.68)]; their active-pole hulls enclose their values by
nonnegativity and partition of unity [PMC09 §1.4.1, (1.58)-(1.61)]. Each hull is intersected with interval
Cox-de Boor evaluation on the current rectangle.

Homogeneous poles, location matrices, knot differences, quotient derivatives and fundamental forms all use
outward-rounded interval arithmetic, with no prior conversion or knot insertion in ordinary floating point
[PMC09 §4.9.1].

Let I be the first fundamental form and II the second form against `S_u x S_v`; principal curvatures are the
generalized eigenvalues of (II, I) [PMC09 §3.4, (3.41)-(3.43)]. **Derived here:** the normal sweep `S + s n`, `0
<= s <= d` (signed), has no focal point on a regular parameter rectangle precisely when `I - d II` is positive
definite, since its eigenvalues in an I-orthonormal frame are `1 - d k_i`.

A negative diagonal entry or determinant proves an eigenvalue has passed zero somewhere along the segment;
Sylvester's criterion proves the positive case without subtracting nearly equal curvature roots. The signed
curvature convention is opposite to the hyperbook's offset chapter; the focal condition is the same [PMC09
§11.3.1, (11.18)-(11.25)]. Face orientation reverses d. A crossing witness must lie strictly inside the trimmed
face; it refuses thickening with kind `Curvature`, naming that source face, before copying or building.

A regular verdict proves only local regularity within the support's knot rectangles: not separation of distinct
normal segments, nor continuity across knots. Periodic or other supports, zero derivative denominators, singular
parameter regions, exact equality and exhausted work remain `Unknown`, and a loose interval never causes a
refusal. The construction runs for unknowns, so complete general self-intersection exclusion is open. The
breadth-first subdivision budgets (4096 rectangles, depth 24) are work choices, not geometric thresholds.

### Elementary focal distances

The offset Jacobian loses rank at a principal focal distance [PMC09 §11.3.1, (11.24)-(11.25); §11.3.2, (11.28)];
with this implementation's normal convention the factor is `1 - d k`. For a sphere or cylinder of radius R the
inward focal distance is exactly R, and these two supports are classified by comparing the signed physical
distance with the stored radius, without the relative allowance used for fitted blends (a radius-10000 sphere
moved inward by 9999.5 remains a radius-0.5 sphere).

For a cone of semi-angle alpha, `r(v) = R + v sin(alpha)` and the principal curvatures are `0` and
`-cos(alpha)/|r(v)|`; endpoint radii and a possible zero crossing give the full focal range over the trimmed v
interval without interior sampling.

For a ring torus (`R > r`) the principal curvatures are `-1/r` and `-cos(v)/(R + r cos(v))`: the inward crossing
is at r everywhere; for outward d the other crossing is `cos(v) <= -R/(r+d)`, and endpoints with the contained
multiples of pi give the cosine extrema over the whole v interval (curvatures with the same normal convention:
[Sage10 principal_directions, torus example]). These reductions are derived here from the parametrizations and
the first and second fundamental forms; arithmetic and trigonometric evaluation keep floating-point limits.

Horn and spindle tori keep strict sampled classification. Classification of general fitted supports stays
sampled with the approximation allowance: not a proven curvature-error enclosure. A general enclosure separating
uncertainty from a strict crossing is open.

### Regularity of network surfaces

A lattice sized from the number of knots does not cover short knot intervals, so rejection witnesses include
every actual knot interval. Acceptance uses the outward-rounded Cox-de Boor derivative enclosure on every knot
rectangle implemented by `focalSweep` at distance zero. The derived identity `det(I) = |S_u x S_v|^2` lets the
cross-product enclosure establish regularity directly without dependency loss in `E*G-F*F`; a zero-free interval
normal proves local immersion and unresolved cells receive no such verdict.

Basis and convex-hull bounds: [PMC09 §1.4.2-1.4.3]. Periodic finite covers use the stored extended flat knot
sequence and cyclic pole copies, matching OCCT's BSplSLib evaluation without fitting. This proves local
immersion only; a global injectivity certificate is open (the argument analyzer can flag valid periodic boundary
seams and does not establish self-intersection freedom within one face).

To reduce interval dependency for rational surfaces, translate homogeneous coordinates by a constant local
anchor, `(X,W) -> (X-aW,W)`: this changes the surface to `S-a` and leaves every first and second derivative
unchanged (derived). All subtraction and multiplication are outward rounded. At distance zero only first
derivatives are needed. Local-immersion queries have no work-count cutoff; the finite subdivision precision
limit reports unresolved, never regular.

### A circular tube's parallel body after its hole closes

For a complete ring torus with core circle C of radius R and tube radius `r < R`, **derived here**: `T = C +
B_r`, so its outward distance-d parallel body is `T + B_d = C + B_(r+d)`, by associativity of the Minkowski sum
and `B_r + B_d = B_(r+d)`. With `a = r + d`, at radial coordinate `rho >= 0` its meridian satisfies `(rho-R)^2 +
z^2 <= a^2`, since the closest point on C has the same azimuth. When `a > R` this is the circle's disk clipped
at the rotation axis.

Revolving the circular arc through `(R+a, 0)` between `(0, +/-sqrt((a-R)(a+R)))`, together with its axis
segment, constructs the boundary without the self-crossing continuation: no surface fitting or new tolerance is
involved, and the exact hole-closure equality `a = R` stays on the engine's route. It applies to a valid
positive-volume solid with one toroidal face, which maps to the new skin; validity, positive volume and the
clearance check still apply.

The volume oracle is the meridian's first moment, `V = 4 pi [ z0^3/3 + R (pi a^2/4 + R z0/2 + a^2 asin(R/a)/2)
]`, `z0 = sqrt((a-R)(a+R))`. This is a trimmed parallel *body*, whereas thickening keeps the input sheet's
normal sweep; trimming arbitrary spline offsets requires tracing all self-intersection components [PMC09
§11.3.2, (11.28)-(11.29)].

### Offsets of moved faces: domains, seams, pipes

#### Extending the approximation domain

The kept-face offset approximates the exact normal offset on the UV rectangle of the constructed face, including
any extension needed to intersect adjacent moved supports. It supplies that rectangle through
`GeomAdaptor_Surface`, whose bounds select an evaluation domain, not `Geom_RectangularTrimmedSurface`, which
rejects bounds outside a bounded basis surface. Approximation tolerance, degree and segment limits are
unchanged.

The normal offset is [PMC09 §11.3.1, (11.17)] and a spline as polynomial spans [PMC09 §1.4.2, (1.62)]. **Derived
here:** evaluating the end span beyond its knot interval continues the same parametric support used by the
extended edge curves; it neither changes the offset definition nor justifies extrapolating a fitted result
outside its certified rectangle. Approximation-error, topology, clearance and self-interference checks still
apply to that rectangle.

#### Both branches of a periodic seam

A seam has two curves on its face but one 3D edge; they may differ by a period in either surface parameter. When
the end faces move, a finite seam pcurve can need the same end-span continuation as a tangent edge. Both
branches are fitted independently at the one-tenth-tolerance interpolation bound, and the modifier gets the
branch selected by that edge use's orientation; binding both curves by face identity alone would overwrite one
with the other. Branches needing no continuation keep their route.

The representation contract is OCCT's `BRep_Builder::UpdateEdge(E,C1,C2,F,Tol)` for a closed face and
`BRepTools_Modifier`'s separate `NewCurve2d` calls for the two uses. The continuation argument is **derived**
from the parameter curves and (11.17). It does not prove that arbitrary continuation is regular or
non-self-intersecting; the constructed body is checked. Closed form for a cylinder with finite spline seam
pcurves: volume `pi (r+d)^2 (h+2d)`.

#### Exact pipes

For an exact pipe `P(s,t) = c(s) + r e(s,t)`, e the unit radial vector in the normal plane of the regular
spine, the normal is radial while the surface is regular; (11.17) gives `c(s) + (r +/- d) e(s,t)`, a pipe about the same
spine (**derived** from the parameterization and normal [PMC09 §11.6.2, (11.111)-(11.115)]). The local condition is that
the shifted radius times the maximum spine curvature is less than one; global separation also matters [PMC09
§11.6.3.4, Thm 11.6.1]. A fitted tensor spline through circular sections is not thereby certified to be this pipe, and
the kept-face operation has no retained spine/radius descriptor for a general spline support; inferring one from a
nominal size and replacing the support would change the input surface and is not done.

### Clearance coverage and rims

Clearance carries separate coverage and resolution states. Empty input, failed extrema, or exclusion of every
region cannot become successful infinite clearance. Distances are measured between whole trimmed regions and the
original skin, boundaries included; an interior lattice alone gives no lower bound. Near openings a region is
measured clear, or wholly inside the opening exclusion, or subdivided.

Work exhaustion is unresolved, never a measurement. *Exclusion bound (derived).* For a closed set A the triangle
inequality gives `d(x,A) <= d(c,A) + |x-c|`; a containing box supplies c and a radius covering the whole trimmed
region, so `d(c,A) + radius < wall` excludes the region. The radius is rounded upwards and the modelling
confusion added. Trimming and distance computations keep OCCT's numerical tolerances: a coverage check using
those solvers, not interval certification of their extrema.

A failed whole-region solve can instead use the lower bound `d(x,skin) >= d(c,skin) - radius`, subdividing until
it resolves the requested clearance; a failed point solve or exhausted subdivision stays unresolved. Distance
queries use independent copies of each trimmed surface.

Spatial edge curves inherited from an import can deviate from the surface within their stored tolerance, and
measuring them would measure representation error as wall thickness; query curves are rebuilt from each face's
pcurves and support with `GeomLib::BuildCurve3d`, whose maximum deviation is returned (unresolved or excessive
deviation fails the query).

For query curve error at most e on either shape, `|d(A,B) - d(A',B')| <= 2e` (derived by the triangle
inequality); subtract it from the numerical distance and add it to an opening exclusion's upper bound, with e =
one tenth of the modelling confusion (a numerical estimate, not a Bernstein-certified enclosure).

A stored edge tolerance [OCCT8 Boolean Operations, "Interferences"] is not a justified allowance to reduce the
requested wall thickness. *Walls and opening rims.* A shell's clearance domain is the inward descendants of the
retained input faces; descendants reached only from removed input faces close the opening and are rims, not
offset walls [OCCT8 BRepOffsetAPI_MakeThickSolid, `MakeThickSolidByJoin`: the first shell holds the input faces
except the removed ones, the second is its offset].

This is tracked by input-face identity through construction, not by surface type, coordinates, a distance
threshold or a hash; a descendant shared with any retained input face is still measured; opening descendants
stay in the history; an empty retained domain still fails coverage.

An opening in a rectangular pocket can give an annular rim slightly beyond the removed face, its outer trim at
the requested offset and its inner trim touching the original pocket walls: independent whole-face distances
give zero there, no positive lower bound exists, and a uniform lattice can miss that zero.

For plane separation h and in-plane distance r from the rectangular opening, `d_open^2 = r^2 + h^2` while the
adjacent trimmed wall can be at r; a spherical exclusion of radius t retains distances down to `sqrt(t^2 -
h^2)`, not t (derived by Pythagoras, `h < t`). Subdivision cannot repair a difference of domains. These are
numerical distance checks [OCCT8 BRepExtrema_DistShapeShape], not interval proofs of extrema.

### Neighbour-bounded planar push-pull

[ZouFeng19 §5.2, Alg. 3] build auxiliary volumes bounded by the neighbours' extended underlying surfaces and the
moved face's start and end supports; a prism is unsuitable for a slanted neighbour. The convex planar
specialization is derived here: write the body as an intersection of outward support inequalities `n_i.x <=
c_i`.

Translating the selected faces by v changes only their constants to `c_i + n_i.v`, and the result is their new
halfspace intersection, exactly the neighbour-bounded volume update; active constraints may appear or disappear
without a topology rule. Moving constants leaves the recession cone unchanged, so a bounded original stays
bounded whenever the target is nonempty.

Triples of independent support planes contain every vertex, giving the temporary box; no artificial box face may
survive, and each surviving face is tracked through Boolean history to its defining support. Both Boolean
differences between the original body and its reconstructed support intersection must contain no solid (this
rejects nonconvex input), and the final topology must pass the validity check. The construction is attempted
only after the sweep rejects a slanted neighbour. Curved supports and nonconvex neighbour-bounded volumes are
open.

### Face-offset endpoint certificate limits

A normal-offset surface [PMC09 §11.3.1, (11.17)] and a valid BREP do not by themselves give a bounded
face-offset result: the boundary must also realize the required incidences within an independently assigned
budget. For every oriented endpoint occurrence measure `|V - C(t_V)|`, and compare the 3D edge curves with both
lifted pcurves. A metadata tolerance enlarged afterwards to cover a measured gap is not an independent result
bound. Closed edges need both endpoint orientations.

The certificate checks every oriented vertex/edge occurrence against `max(operation tolerance, input vertex
tolerance, carried input vertex tolerance + incident normal spread)`, assigned before any incident curve is
fitted; both lifted seam pcurves are checked after rebuilding and optional surface approximation. The inherited
normal-spread estimate is sampled, not an interval enclosure.

In the fixed-incidence construction, continuing a pcurve beyond an approximately tangent seam's old trim need not keep
it on the other offset support, and a local point-on-support solve may locate a vertex while its incident curves still
miss it (an SVD Gauss-Newton vertex refinement returned valid bodies with gaps 700 to 1000 times the input tolerances).
Two representation defects were found:
1. *Seam carried through a moved junction.* An inward junction on a cylindrical support moves in azimuth while its chart
   seam stays fixed. A seam is not a physical intersection: carry both pcurves through the projected endpoint UVs with the
   same linear endpoint correction, then fit their common lift. **Derived:** if the original pcurves differ by a period
   vector kP, adding the same correction `delta(t)` preserves `S(q(t)+delta(t)) = S(q(t)+kP+delta(t))` on a periodic
   carrier; the support and physical trims stay fixed. This needs regular endpoint projections and a valid cut in the
   retained region (validity, normal distance and interference stay separate checks) and is not a general existence
   proof for arbitrary seam routing.
2. *Closed edge without periodic parameter.* A closed topological edge need not have a periodic curve parameter
   [OCCT8 BRep_Tool `IsClosed`, `CurveOnSurface`: closure is by shared vertices; the two seam pcurves are distinct].
   Reusing its old parameter width for a moved endpoint can put the final section off the junction; where the candidate
   final section exceeds the preassigned endpoint budget, the final occurrence's section parameter is solved
   independently. A defect with a seam as retained cut (normal spread exactly zero) is held to the operation's 1e-6
   allowance and is not an unresolved angular-budget case.

#### Why an old trimmed-contact bound does not certify a continued contact

Source lifts continued independently (by homogeneous de Casteljau continuation of the end profile span, then
rotation, which matches the offset `S + d n` to floating-point error) already differ at the first refused
contact by 3e-5 with normals differing by 1.6e-3 although they agree at the old vertex. The triangle inequality
gives `gap_after <= gap_before + |d| |n1-n2|` (derived). It neither authorizes a larger output tolerance nor
shows that a jointly constructed continuation cannot exist.

Closest-point stationarity: [PMC09 §7.2.1, (7.8)-(7.9)]; SVD least squares and Gauss-Newton: [NocedalWright06
§10.2-10.3, (10.20)-(10.26)].

Neither establishes that a locally refined vertex connects to the desired intersection branch, a failed
iteration is no impossibility proof, and reconstructing nearby exact section intersections and minimizing
contact distances in model space also failed the measured continuation; none shows that no offset exists within
the input uncertainty. *Structural alternative (proposed, not implemented).* Construct each junction together
with its incident contact curves, track intersection branches before trimming, and carry explicit endpoint and
edge/support budgets through the fit. [Mukundan04 §2.2, (4)-(6), §3-4, Thm 1] distinguish isolated tangential
contacts, branching and higher-order contacts, then validate tracing and map parameter enclosures to model-space
enclosures; their assumptions (an identified branch, a bounded starting point, regular rational interval
surfaces) make this no ready theorem for arbitrary normal-offset surfaces or approximate contact seams.

The [OCCT8 BRepOffsetAPI_MakeOffsetShape `PerformByJoin`] contract distinguishes intersection joins from added
arc faces and promises no self-intersection elimination. A general construction must establish its own existence
conditions and final bounds.

### Separated tangent junctions

The fixed-incidence construction runs first and keeps its successful results. After a vertex fails the
point-on-every-support check, a single-shell solid with three incident faces and a retained tangent edge can
instead split the junction. The tangent edge ends at its original endpoint transported to the offset; each of
its two continuations intersects the third offset support independently, and the two sharp edges use those
separate endpoints. The old support surfaces stay exact.

A ruled strip `R(u,v) = (1-v) C0(u) + v C1(u)`, boundary curves on [0,1] [PMC09 §9.7.1, (9.74)-(9.76)], fills
the boundary between the continuations and its transverse end closes the third face; face lineage gives the
strip all three incident parents. `BRepFill::Face` builds it through `GeomFill_Generator`.

It is a boundary closure, not a claim that the strip is a regular developable surface or an exact
normal-distance envelope (normal offset: [PMC09 §11.3.1, (11.17)]; separately generated edge and corner patches
in solid offsets are discussed in §11.1.4, which does not prove this strip's validity). Unsupported arrangements
refuse. *Checks.* The seam's two transported endpoints must agree within the operation tolerance; rail fits use
one tenth of it; the transverse chord must project to the third support within it.

Every oriented endpoint and each new boundary's lifted pcurve is checked against that bound independently of
output metadata, and validity, self-interference, signed-distance, volume and whole-body clearance checks still
run.

Curve/support checks are adaptive fits followed by sampled validation, not interval certificates of the whole
strip. *Shell orientation.* Shell reconstruction keeps the orientation of each source-face occurrence, including
faces outside the split junction: the modifier's shape map identifies shapes independently of occurrence
orientation, and its stored face is not a replacement for the oriented use [OCCT8 Modeling Data, "Orientation";
`BRepTools_Modifier::Rebuild`, "copy the sub-elements"].

The stored orientation of an unchanged face can give a shared edge two identical uses, so a closed collection of
individually valid faces is not an oriented boundary. Carriers, curve parameters and contact bounds are
unchanged. **Derived diagnostic:** for two regular contact continuations with a common position and tangent, the
normal separation is `0.5 (k0-k1) s^2 + o(s^2)` in the common tangent arc coordinate s.

The curvature is taken on the continued side of a boundary knot (a symmetric finite difference across a
second-derivative break averages the two sides). This Taylor estimate diagnoses the failed old incidence and
never assigns a result tolerance; it does not certify existence, regularity or nonintersection of a finite
strip.

### Transported contact uncertainty is confined to its measured domain

Let `p(t), q(t)` be corresponding lifts of a source edge and `n(t), m(t)` their outward unit normals. From the normal-offset
definition [PMC09 §11.3.1, (11.17)] the triangle inequality gives the **derived** bound for equal signed moves d:

`|(p+d n) - (q+d m)| <= |p-q| + |d| |n-m|`.

For unequal signed moves replace the second term by `|d1 n - d2 m|`; if the representative 3D curve rather than
the other lift is used, its discrepancy is added. A bound on the original closed parameter interval does not
bound an extrapolation outside it, so a large continued residual cannot justify itself by replacing the source
uncertainty; and a measured nonzero input error is no necessary lower limit on output accuracy (the inequality
is an upper bound, not a nonexistence theorem).

The operation's 1e-6 allowance is not replaced by a constant near 1e-5: on 47 generated sources with approximate
tangent contact against revolution supports (input edge tolerances 1e-7) the sample maxima were 5.2e-8 for the
source lift disagreement and 5.5e-8 for the pointwise transported bound at `d = 0.05`; sample maxima are
diagnostics, not certified suprema, and the endpoint and paired-pcurve checks remain necessary.

A reference result preserving the whole face-edge-vertex incidence graph, or returned without error, is no
strict geometric oracle: bodies were seen with curve endpoints within 7.6e-10 whose junction points were 2.3e-5
from the exact offset revolution support, and others with endpoint discrepancies above 5e-5. *Offset of a
revolution by its meridian (derived).* Write the generating curve as `(r(v), h(v))` about the axis, `r > 0`, `L
= sqrt(r'^2 + h'^2) > 0`.

For surface-normal offset d the offset meridian is `(r + d h'/L, h - d r'/L)`, by evaluating `S_u x S_v` for a
revolution. For query cylindrical coordinates `(rho, z)`, minimize the Euclidean distance to that curve; angular
minimization aligns the radial directions while the offset radius stays positive. Scalar minimization over an
enlarged parameter interval agrees with direct 3D projection to about 2e-15 on the first rejected junction: a
cross-check of that measurement, no whole-domain existence or global-minimum certificate.

Exact versus approximated offsets: [PMC09 §11.3.7]. A continuous contact construction still needs branch
ownership and endpoint compatibility; the initial-component assumptions of [Mukundan04 §1, 2.2, 4] are not
supplied by an independently continued approximate contact, and neither the reference observation nor the
finite-point bound closes that problem.
