<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepBlendExt and BRepFilletExt: rolling-ball blends

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Problem](#problem)
- [Construction](#construction)
- [Convexity](#convexity)
- [Sources read](#sources-read)
- [Connected centre-length fitting](#connected-centre-length-fitting)
- [Support knots of continuity C2](#support-knots-of-continuity-c2)
- [Arc length beyond a trimmed rail](#arc-length-beyond-a-trimmed-rail)
- [Fixed-contact conic windows](#fixed-contact-conic-windows)
- [Rigid-motion invariant collision prefilter](#rigid-motion-invariant-collision-prefilter)
- [Separation certificate before the exact collision distance](#separation-certificate-before-the-exact-collision-distance)
- [Contacts beyond an entire smooth support](#contacts-beyond-an-entire-smooth-support)
- [Closed offset-intersection station carriers](#closed-offset-intersection-station-carriers)
- [Open offset spines walked by arc length](#open-offset-spines-walked-by-arc-length)
- [When the spine is used](#when-the-spine-is-used)
- [Radius-law extrema](#radius-law-extrema)
- [Steep and flat wedges](#steep-and-flat-wedges)
- [A size that changes along the edge](#a-size-that-changes-along-the-edge)
- [Why spatial stations stay](#why-spatial-stations-stay)
- [Connected limited contours](#connected-limited-contours)
- [Over a side narrower than the blend, and a pinched body](#over-a-side-narrower-than-the-blend-and-a-pinched-body)
- [A run back of no length](#a-run-back-of-no-length)
- [Where a contact folds on any face](#where-a-contact-folds-on-any-face)
- [Where a round crosses itself between its contacts](#where-a-round-crosses-itself-between-its-contacts)
- [Further constructions](#further-constructions)
- [Unproven choices](#unproven-choices)
- [Stock counterpart](#stock-counterpart)
- [Local refinement without an iteration budget](#local-refinement-without-an-iteration-budget)
- [Located continuity events in a fitted blend](#located-continuity-events-in-a-fitted-blend)
- [Construction boundary roles](#construction-boundary-roles)
- [Compatibility before a second-order join](#compatibility-before-a-second-order-join)
- [Derived negative control](#derived-negative-control)
- [A remaining degree of freedom in symmetric quintic sections](#a-remaining-degree-of-freedom-in-symmetric-quintic-sections)
- [Internal curvature joins and boundary-ribbon compatibility](#internal-curvature-joins-and-boundary-ribbon-compatibility)
- [Contact ownership on a joined blend sheet](#contact-ownership-on-a-joined-blend-sheet)
- [Centre-path arc coordinates for a defective sheet](#centre-path-arc-coordinates-for-a-defective-sheet)
- [Source-normal events of a centre carrier](#source-normal-events-of-a-centre-carrier)
- [Continuity intent after limiting a blend](#continuity-intent-after-limiting-a-blend)
- [Smooth section-policy selections](#smooth-section-policy-selections)
- [Experimental terminal restriction continuation](#experimental-terminal-restriction-continuation)
- [Boundary intervals after two support crossings](#boundary-intervals-after-two-support-crossings)
- [Non-isoparametric planar stops across spline breaks](#non-isoparametric-planar-stops-across-spline-breaks)
- [Carried-sheet corner boundaries](#carried-sheet-corner-boundaries)
- [Stops reached past the edge's vertex](#stops-reached-past-the-edges-vertex)
- [Split corners of more than three conic blends](#split-corners-of-more-than-three-conic-blends)

## Notes by topic

### Problem

Two supports S1, S2 meet along an edge; for a radius r (constant, or a law r(s) along the edge) find the
envelope of balls of radius r that touch both. The ball centre c lies on both offset surfaces, c in O1(r) ∩
O2(r) with O_i = S_i + r n_i (the sign follows convexity); the contact points are p_i = c - r n_i(p_i), and the
blend is swept by a section curve between p1 and p2 in the plane normal to the centre curve. G1 holds at p_i by
construction (the section is tangent to the ball).

The offset surface is regular where 1 - r k_i is non-zero for the principal curvatures k_i of the support on the
side the ball is on, so a radius reaching the smallest concave radius of curvature is refused rather than built
[Jaiswal et al.].

A section curvature-continuous with its support needs only the normal curvature of the section to match that of
the support in one direction across the join, by the Linkage Curve Theorem [Salvi and Várady], so the section's
third pole is set by the support's curvature; that rule is a section *policy*, which a caller supplies (see
[Section policy](../../../README.md#section-policy)).

### Construction

The centre curve is marched from a start point by Newton correction of the two
offset equations with curvature-controlled steps [Jaiswal et al., after Patrikalakis and Maekawa];
exact sections are fitted by a spline with knots at the data's continuity breaks and interpolation
at Greville abscissae (Schoenberg-Whitney) [de Boor]; the sheet is attached by trimming the supports
to the contact curves and the stopping faces to the ends. The conic section is the rational quadratic
Bézier whose shoulder weight is sin(w/2) rho/(1 - rho) for a wedge angle w, which is the circular arc at
rho = 1/2 [standard conic form]. Supports are extended past their ends by continuing the last
polynomial span, which keeps shape and parametrisation [Hu, as restated by Zang et al.].

The source follows these stages, one file each, with each stage's steps as named functions:

| File | Stage |
|---|---|
| `BRepBlendExt_Solver.cxx` | The request: the routes that lift features off, replay overruns or blend one kind of edge first |
| `BRepBlendExt_Stages.cxx` | The pieces, their joints, their sheets (`Blend::sheets`) and the splits where a contact crosses an edge |
| `BRepBlendExt_Attach.cxx` | Attachment: the chains' ends and the assembly |
| `BRepBlendExt_AttachStops.cxx` | The faces that stop a blend (`Attach::stop`, `Attach::stops`) and notches |
| `BRepBlendExt_AttachMiters.cxx` | Mitres, and balls rolled round a re-entrant vertex |
| `BRepBlendExt_AttachCorners.cxx` | Corners of three or more blends (`Attach::CentralSplit`) |
| `BRepBlendExt_AttachFolds.cxx` | Loops where a contact folds over itself, and their patches |
| `BRepBlendExt_AttachFaces.cxx` | Contacts, blend faces, continuity, tolerances and collisions |

### Convexity

The sign of the offsets is the sign of the turn h = (n1 x n2) . d/|d| = +-sin(theta) along the edge, proven on
every part of it (BRepBlendExt_Convexity.cxx).

With the unnormalised normals X = S_u x S_v and the tangent d = d/dt S(u(t), v(t)) as polynomials in one
parameter, H = (X1 x X2) . d is a polynomial with h = H/(|X1| |X2| |d|); its Bernstein coefficients bound it
with the correlation between the faces kept (faces tangent along the edge give coefficients at rounding), |X|^2
is bounded below by its least coefficient, and halving a part brings the coefficients to the values
quadratically [Lane and Riesenfeld 1981].

Rational pieces compose homogeneously; a quadric's or torus's normal is the gradient of its implicit equation at
the edge's point (a point within the edge's tolerance of it), a plane's is constant. A part is halved until it
is proven of the edge's sign or tangent all over (|h| below the sine of the tangency angle for the blend's
size), or its middle shows the other sign: an edge with convex and concave parts is refused by name, since the
ball would change sides where the faces are tangent between them.

Certified for spline and Bezier surfaces and extrusions of spline curves (with spline, Bezier or line parameter
curves), planes, cylinders, cones, spheres and tori; revolutions and offsets are read at the places k/16 of the
edge.

Where the faces are tangent at some places and turn one way at the others, the blend narrows to the edge's point
there: the march starts where the faces turn most if they are tangent at the middle, a station where the normals
agree to the angular resolution is the ball resting on both at the edge's point (as at a tangent end), a section
whose tangents at the feet are opposite to within sqrt(1e-12) takes the middle of the feet as its apex (the
limit as the wedge opens to a half turn; derived), and the ball's Newton solve accepts a residual that stops
falling below its forward error bound: with the faces' first derivatives off by eta <= gamma_k (k their
evaluation's operation count, Higham 2002, 3.1), the unit normals are off by 2 (2 eta + gamma_2) |S_u||S_v|/|S_u
x S_v| each, and the plane residual arm . (n1 x n2)/|n1 x n2| by |arm| 2 e/|n1 x n2| + gamma_3 (|arm| + |c| + |E|)
with e the two normals' errors and sqrt(6) gamma_2 (derived; no chosen constant), and only where the
section's width times the plane's uncertainty 2e/|n1 x n2| stays within the solve's tolerance (a narrow section
between all but tangent faces; not a ball held between faces turned against each other).

### Sources read

Jaiswal, Ghosal, Gurumoorthy, "Circular blending", DETC 2001 (theory, algorithm,
results; Hermann 1995 and Lukacs 1998 are cited second-hand there); Salvi and Várady, "Multi-sided
surfaces with curvature continuity" (2014; statement of the Linkage Curve Theorem, appendix proof not
read); Zang, Liu, Lai, GMP 2008, sections 1-3 (Hu's extension; Hu 2002 not obtained); de Boor 1972.

### Connected centre-length fitting

At a defective rail knot in a constant-size
curvature-continuous blend, connected pieces now share one numerical centre
carrier. Each piece uses its restriction of that carrier's integrated length;
the joined parameter widths are proportional to those lengths. Transforming a
retained terminal jet uses the chain rule, including the second derivative of
the inverse length map. Split-normal coupling uses the same derivative rates.
Independently normalising the pieces without these rate changes does not
preserve their coupled tangent planes.

The centre carrier retains a C0 knot at each change of rolling support: the
pieces share a centre there, but equal centre tangents have not been established.
This is a spline-space requirement, derived from the continuity of the data;
smooth refinement alone does not represent that event. Carrier refinement is
not refused at a fixed coefficient count. It retains the station-solve length
tolerance and requires new representable knots until the sampled residual fits.
Coefficient count is a representation cost, not a geometric obstruction.
The independent two-line centre family checks position, exact piecewise arc
length and its inverse, under rigid transformations and scales 1e-3 through 1e3.

Where the section law has no framed break, the two terminal homogeneous jets
are averaged in physical centre length. Quintic Hermite displacements carry
these first and second jets, vanishing through order two at their inner window
boundaries. The existing radius-bounded coupling windows cannot overlap in
the interior of a piece. Rational weights are included, so a homogeneous scale
change does not change the construction. This jet transport and the polynomial
window formulas are derived applications of the chain rule and Hermite
interpolation. They do not prove that an arbitrary prescribed section law is
C2 across a support event; such events retain their knots and framed law.

### Support knots of continuity C2

A curvature-continuous section follows the support's curvature, so where a foot crosses a knot at which a spline
support is only C2 the exact law is C0 in its poles: the sheet's second form jumps there.

Within a window of half-width h in the support's coordinate x from the knot, the section is that of both
adjacent polynomial spans, each continued with its own station solve and frame, as P_a + S(z)(P_b - P_a) in
homogeneous coordinates, S(z) = 6z^5 - 15z^4 + 10z^3 and z = (x + h)/2h [Perlin, *Implementing Improved Perlin
Noise*, GPU Gems (2004), section 5.4]. S has zero first and second derivatives at both ends, so the law is C2
there and the fit takes knots of multiplicity three at the window's ends and at the knot.

The window is derived, not tuned (derived here): the spans agree through order two at the knot, so on either
half the law departs from its own exact branch by S(z)|D(x)|, D = P_b - P_a. Where |D(x)| <= |x| m, the
departure is at most c m h with c = max r S((1-r)/2) = 0.0705529 (r solving 4 - 15r + 20r^3 - 9r^5 = 0); this
bounds the poles, hence the section, and the section's curvature at the foot, which is the contact's across the
blend.

The feet part as |x|^3 m', and the contact's normal curvature along the blend departs by at most k m' h / v^2, v
the feet's speed in x, k = max |(S(z)(2z-1)^3)''|/4 = 0.328189 over [0, 1/2]. h is the widest window (within
half of each adjacent span) whose departures, read from both spans at x = -h and h, keep the section within its
limit less the fit tolerance and the contact within the fit's own limit on its curvature (limit::bend). These
are leading-order bounds: m and m' are secants over the window, not certified suprema.

A local split join (contacts whose supports agree through curvature) whose parts still bend apart by more than a
delivered seam receives the connected centre-length fitting above, with second jets coupled in the section's
interior only: the first and last three homogeneous columns of a clamped section fix both contact jets through
order two, so the contacts stay exactly as fitted.

Sources: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided
Design and Manufacturing* (2009), section 2.1, equations 2.6--2.11; Peters,
*Geometric Continuity* (2001), sections 2--3. Centre fitting and geometric
acceptance remain sampled numerical checks. The generated polynomial and
rational-sheet family checks common terminal jets through order two, unchanged
window exteriors, rigid transforms and scales 1e-3 through 1e3.

The internal jet-transport primitive also provides contact-preserving coupling
for a clamped transverse basis. Its first and last three homogeneous coefficient
columns determine the value, first derivative and second derivative at each
contact. Leaving these columns unchanged therefore preserves the complete
contact jets throughout the longitudinal window, including mixed derivatives
and rational weights. Only interior columns receive the split-jet correction.
This is a derived endpoint identity, not a proof of compatibility in the
preserved ribbons: their remaining split residual still needs checking. The
ordinary production caller retains its existing coupling. The generated family
also checks all six components of the contact jets before and after transport.

### Arc length beyond a trimmed rail

A rail can outlive the edge that supplied its parameter
range: the adjoining supports still carry the ball beyond an endpoint. For a nonperiodic multi-span
curve, `GCPnts_AbscissaPoint::Length` sums only the intersections with the adaptor's continuity
intervals (`computeLength`, composite case); an interval wholly beyond the edge consequently has
zero reported length. Expand a copy of the adaptor's parameter range before asking for arc
length: its original spline and rigid transform are retained, including knots hidden by the edge
trim. Integrating one supposed end span without exposing those hidden knots is insufficient.
The curve itself is neither fitted nor extended by a new geometric construction.

For the inverse, use the unchanged solver inside the edge's range. Outside, start with the local speed estimate
`distance / |C'(t)|` and double the parameter interval until its integrated length brackets the requested
distance. This bracketing argument is derived here: on a regular branch, `d(length)/dt = |C'(t)| > 0`, so the
length is strictly monotone; if a finite parameter realizes the requested length, a finite number of doublings
reaches it.

Non-finite arithmetic and a step with no representable parameter progress do not supply a bracket and are
refused. Outside the trim, integration and inversion use the existing scale-resolved `limit::solved` tolerance;
no seam or construction bound changes. Integration and root solves remain numerical, without an interval
certificate; no global existence claim is made past a singularity or where the continued curve has insufficient
length.

The closed-form check uses (t,t^2), whose primitive is
`t sqrt(1+4t^2)/2 + asinh(2t)/4`, across both endpoints and an inserted knot, both directions,
rotated and at scales 1e-3, 1 and 1e3 (294 interval/inverse checks). This primitive is derived here.
A held-out rational quadratic circular arc, also with an inserted knot and trimmed endpoints,
uses radius times the difference of endpoint `atan2` angles as its oracle (147 more interval/inverse
checks at those scales). A C1 curve that is a line up to t=1/4 and (t,(t-1/4)^2) beyond it tests
both trims that hide that knot, in both directions (98 interval/inverse checks); the earlier
single-span extension integration fails this family. No oracle uses the integration code for its
expected length.

This is the arc length of the **rail**, not a replacement of the centre path by an arc-length march.
Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing*, hyperbook
edition (2009), [section 5.8.2.3, equations
(5.101)-(5.105)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html), give the normalized
intersection tangent and parameter ODEs; transverse intersection is required, and the normalized cross product
fails at tangency.

Their topology discussion also requires a start on every connected branch before marching. Jaiswal, Ghosal and
Gurumoorthy, [*An Efficient Algorithm for Circular Blending of
Edges*](https://mecheng.iisc.ac.in/~asitava/sharad.pdf) (2001), introduction and theoretical background,
explicitly separate geometric spine construction from trimming and integration into a valid body: that paper
addresses the former.

Neither result by itself proves ownership of every retained trimmed region or that a coupled section satisfies
its prescribed law; those remain separate acceptance requirements.

### Fixed-contact conic windows

A final closed offset-circuit route can now attach a non-circular conic skin for a
cylindrical cavity meeting a transverse plane, at constant circular radius
`0 < r < R`. It follows every previous whole-request route. Neither previously
successful requests nor the crossing and solid acceptance checks are bypassed.
The carrier parameter frame is immutable: restricting a piece at a support seam
must restrict that circuit, rather than map a new complete circuit onto the piece.

Let `q` be the cylindrical contact's signed phase from the acute radial direction, `A=R+r`, `h=2 asin(r/A)`, and
`x=q/h`. The complete ellipse double-point domain has half-width `asin(sqrt(((r/A)^2-c^2)/(1-c^2)))`, where `c`
is the absolute cosine between the plane normal and cylinder axis. It is enclosed by `h` when `r/A >= c`. These
formulas and the enclosure are derived; they are not a published construction theorem. Outside `|x|<1`, the
circular generator remains unchanged.

Inside, replace its shoulder `M` by `M + (1-x*x)^2 (M0-M)`, with `M0` the intersection of the cylinder-axis line
through its cylindrical contact with the plane. Both shoulders lie in both contact tangent planes. Positions and
tangent planes at both contacts remain exact analytically. The window and its first derivative vanish at the two
ends, so the retained/replacement first jets match.

The existing homogeneous fitter receives these sections and both contact pcurves together; knots of multiplicity
four in its degree-five longitudinal basis represent the actual C1 window breaks. Homogeneous subdivision into
two transverse spans follows the shoulder change.

The ordinary fitter controls error numerically at its refinement sites. Neither
continuous point/derivative error of this spline approximation nor its global
embedding has a runtime interval certificate. The analytic first-jet identities
do not supply those certificates. Delivered results use the same validity,
closed-solid, sampled join and retained-section checks as the other blend routes;
these checks are numerical, not a proof of global injectivity.
The generated cavity/plane tests retain their original rolling-ball section-error
bound outside the replaced sector, and independently test the immutable carrier
under seam restrictions. No acceptance bound changes.

Sources read: Jaiswal, Ghosal and Gurumoorthy, *An Efficient Algorithm for Circular
Blending of Edges* (2001), introduction, pp. 1-2,
<https://mecheng.iisc.ac.in/~asitava/sharad.pdf>, separates geometric spine
construction from trimming and solid integration. Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
§5.8.2.3, equations 5.101-5.105,
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html>, gives
transverse SSI continuation and its branch-start/tangency conditions. These
sources motivate the station/attachment separation; the compact shoulder
homotopy and its identities above are derived here.

### Rigid-motion invariant collision prefilter

The collision stage keeps its world-axis bounding boxes and exact face-distance
query. Where a box pair overlaps, it additionally tests a positive-weight tensor
spline against the other face's entire support plane. The signed plane distances
of all control poles enclose every point of the spline by the convex-hull
property. If this interval, enlarged by coordinate roundoff, lies strictly beyond
the sum of the two face tolerances and both boundary tolerance envelopes on
either side of zero, no point of the two faces can collide. A trimmed face is a subset of that spline, so trimming cannot
invalidate the separation. Other support types and inconclusive intervals retain
the original exact query. No geometry or collision tolerance changes.

The enclosure is the standard NURBS convex-hull property: Alexander Efremov,
*Efficient Ray Tracing of Trimmed NURBS Surfaces*, master's thesis (2004),
§2.5.3, property 3, pp. 23-24,
<https://dcgi.fel.cvut.cz/~havran/ARTICLES/EfremovMasterThesis.pdf>.
The separation inference above is derived. It avoids expensive distance queries
caused solely by rotation of long thin surfaces' world-axis boxes.

### Separation certificate before the exact collision distance

The collision stage asks one yes/no question per pair of a new blend face and a face it is
not next to: do they come within the sum of their tolerances? The exact query
(`BRepExtrema_DistShapeShape`) answers it with a global minimum distance, through a Lipschitz
optimiser on every edge pair. Face boxes cannot help when one face wraps around the other, as
the far half of a closed tube does, so most of the stage's time went to pairs many millimetres
apart. `BRepBlendExt_Separation.cxx` decides those pairs first, and the exact query decides the
rest, so no verdict can change.

- **Cover.** A face lies on its surface inside the parameter box of its edges' curves on it
  (curves bounded by their poles or exactly: lines, circles, B-splines, Beziers). That patch is
  converted exactly into a rational B-spline with positive weights: planes, quadrics, tori,
  B-splines, Beziers, and surfaces of revolution and extrusion of lines, circles, ellipses,
  B-splines and Beziers. An offset surface's points lie within its offset of its basis, so the
  basis is converted and the offset added to the face's growth. Edges and vertices lie within
  the face's largest tolerance of the patch (derived). Bounded directions must hold the box;
  an extrapolated patch is not covered.
- **Hulls.** A positive-weight rational patch lies in the convex hull of its poles
  [PT97 P4.4, tensor-product form]. Two hulls whose pole projections on a unit direction are
  further apart than the tolerance plus both growths and coordinate rounding hold faces that
  are (derived: projection does not lengthen). Tried directions: the coordinate axes (box
  distance), the line between the cells' middle points, and each cell's normal there. The
  normal direction makes the bound converge to second order for nearly parallel pieces.
- **Halving.** Pairs not separated halve the cell with the larger pole box, across its longer
  control polygon, from the whole cover each time (so the rounding is that of one
  segmentation) [JC98].
- **Stops, all derived.** Not proven, and left to the exact query, when two cover points are
  found within the tolerance (the middle points, or a local closest pair from
  `Extrema_GenLocateExtSS` started there when the hulls overlap): no hull of these covers can
  then separate them. Also when both cells' pole boxes are within the tolerance across: their
  points are within three tolerances, which no finer cell resolves. There is no work budget.

Measured on a matrix of blend supports and a set of blend requests (1,505 pairs that reach the exact query on the
unchanged code): 1,474 certified, none of them within the tolerance by the exact query (smallest
certified distance 0.0067 mm); exact query time 84.9 s before, 0.57 s on the 31 pairs left; the
certificate 0.19 s in all. Every result is identical. CPU per request, warm medians:

| Request | Before (ms) | After (ms) | Ratio |
|---|---:|---:|---:|
| Torus beside a B-spline, concave groove, round, r 1 | 1,093 | 222 | 0.20 |
| Same, conic (rho 0.3) | 1,112 | 212 | 0.19 |
| Same groove, curvature-continuous, r 0.2 | 1,825 | 175 | 0.10 |

[PT97] L. Piegl, W. Tiller, *The NURBS Book*, 2nd ed., Springer 1997, P4.4.
[JC98] D. Johnson, E. Cohen, "A framework for efficient minimum distance computations",
ICRA 1998.

### Contacts beyond an entire smooth support

A contact can lie beyond a smooth support boundary for the whole retained rail
interval. There is no longitudinal crossing at which to split that interval.
The continuation retains the source station frame, projects the middle contact
onto the neighbouring support, and solves the same offset equations there. A
projection is only a Newton seed: the solved foot must classify IN or ON against
the finite destination face before that face can become a rolling support. The
existing fit, contact, section and solid checks still apply. A support bounded
only by the source rail and crossed boundary is consumed together with its
boundary edges; adjacent faces are rebuilt on the continued contacts.

After end cuts, mitres and ridge replacements restrict each continued contact to its retained intervals,
acceptance checks the support pcurve attachment actually stores: joined where pieces share a contact, and
projected from the exact sheet where it replaces the fitted endpoint map. Acceptance intersects this UV curve
with every trim boundary and classifies each interval between events, including the contact endpoints. Ends of
tangential overlap segments are events too.

On periodic carriers, conservative UV boxes determine every translated boundary that can intersect the contact;
each witness belongs if it is in the face interior or within the face's own model-space tolerance of a trim
evaluated on the carrier. Boundary distances use point/curve extrema and both curve endpoints, including
degenerate pole edges. The lifted trim is essential: the independently stored 3D edge can differ from it by the
edge tolerance.

This union is the definition of membership in the closed face's tolerance band (derived); it is invariant under
periodic parameter translations and does not turn a length tolerance into an angular tolerance using a speed
measured away from the witness. Membership is constant between exact boundary events outside the tolerance band
(derived from continuity). The event locations and classification remain numerical OCCT operations, rather than
an interval certificate.

Neither a carrier extension nor a missed or incomplete intersection computation establishes finite support
membership. A contact outside its finite support is refused before attachment.

This construction handles whole-interval crossings in the finite-support routing.
Boundary normal agreement is sampled using the existing tangent policy;
it is not a continuous smoothness certificate. The regression family varies
radius, scale and rigid pose and independently measures the ball-contact defect.
Open intrinsic center-path correspondence and overlapping sheets remain separate
work; this construction does not trim them.

For a curvature-continuous request, the retained internal knots also pass the
existing absolute normal and second-form bounds in `sheetJoins`. The source
parameter guide is a coordinate after transfer, not a retained rolling support:
its section-field break cannot excuse an internal curvature jump. Nor does a
small positional sag certify curvature continuity. These continued sheets use
the absolute bounds without the section-law or sag exemptions. One-sided spline
jets are exact on their chosen spans; the seventeen transverse sites per span
are a sampled check, not a continuous certificate. A failed check refuses the
construction; it is not evidence that no compatible curvature-continuous surface
exists. Constructing compatible section second derivatives remains open.

### Closed offset-intersection station carriers

After the existing whole-request constructions refuse, a separate carrier may be obtained by intersecting the
constant-distance offset supports with `GeomAPI_IntSS`. Canonical offsets are represented exactly through
`Geom_OffsetSurface::Surface`; other offsets retain their surface evaluator. The source-edge adaptor and its
parameter range are kept separately. The branch must be a closed circuit through the existing middle station,
with compatible end tangents. Its orientation follows that station.

A closed clamped SSI spline is evaluated by wrapping its finite domain, never by continuing its end polynomial.
Every ball is still corrected on the original offset equations. These numerical tests identify a branch, not a
proof that all branches were found. Open chains, carried ends, variable distances and cliff contacts are not
handled by this carrier construction.

The intersection tangent construction requires transverse normals: Patrikalakis, Maekawa and
Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
[section 5.8.2.3, equations (5.101)-(5.105)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html).
The same section distinguishes tracing from finding a seed on every component. The choice of
one circuit through the middle station and its affine parameter correspondence are implementation
choices, not a published completeness theorem. Independent plane/cylinder equations and winding
checks cover 24 combinations of angle, radius, scale and pose in the native test.

The carrier alone moves folded station maps to the unchanged interior-crossing check; it
does not repair overlapping envelopes. On the cavity/plane class the conic-window construction
above supplies that replacement. The immutable carrier frame keeps restrictions at support seams
from stretching another circuit across a fragment. Forcing a new carrier on completed circuits
can relocate the support seam, so this construction stays after the previous whole-request
routes. Other crossing classes still require their own compatible contact arrangement.

### Open offset spines walked by arc length

The same final route now gives an open constant-distance rail (no size law, no cliff or trimmed side) its
offset-intersection spine.

The spine is traced through the balls themselves: the ball equations c0(u0, v0) = c1(u1, v1), each centre c_i =
S_i + r_i n_i off its support, are three equations in the four foot parameters, and where their Jacobian J (3 x 4)
has full rank their solutions are a curve whose tangent is J's null vector, the signed 3 x 3 minors of J (a 4
x 4 determinant with a repeated row vanishes; **derived**), with the centre moving by J_u0 t_0 + J_v0 t_1.

It is followed by pseudo-arc-length continuation (Allgower and Georg, *Introduction to Numerical Continuation
Methods*, SIAM 2003, chapters 2 and 6): an Euler step along the tangent, scaled to unit centre length, then
Newton on the ball equations within the hyperplane across the tangent through the step's end.

The carrier is the C1 cubic through the traced balls with their unit tangents scaled by the chords between them
(Hermite); a step is kept where that cubic's middle lies within the joining tolerance of the path's own ball
there, found the same way, and the next step is scaled by the fourth root of the error's ratio to the tolerance,
at most doubled, as the cubic's error is of fourth order (a step policy: the balls are corrected on both
supports at every station, so the carrier only places them).

On a span the carrier's squared speed |dc/dtau|^2 is a quartic in the span's own parameter tau, with dc/dtau = A +
B tau + C tau^2 from the span's Bezier poles (**derived**); its arc length is the square root of that quartic
integrated by 8-point Gauss-Legendre (Abramowitz and Stegun, 1964, 25.4; exact for polynomials up to degree 15,
and the speed varies little across a span: the order is a choice, not a proven error bound) and inverted by
Newton within the span.

The general centre-length map (adaptive quadrature over the whole spline) took some 19 us an inversion and some
11,000 inversions on a refused chain, two thirds of that stage's time.

The trace starts at the rail's middle ball and goes each way to the ball at the edge's end there, found exactly
as the march finds its first ball (grown from the edge to its size, through the edge's points or, for a rail
found through its centres, through the centres across it). Those end balls are the correspondence of the spine's
ends to the edge's.

The path reaches an end ball where it crosses the plane through its centre across its path at that centre (the
crossing solved on the path by Newton with that plane as the fourth equation); it crosses there along the ball's
own path direction, so a step, which turns the tangent by less than a right angle, cannot pass it unseen.

It proves the end ball on another branch, so that the spine does not carry the rail, where the path closes
(crosses the middle ball's plane at the middle ball from behind) or leaves the supports' domains first (a
parameter neither periodic nor unbounded past its surface's bound, or, where a support runs on without bound, a
centre farther from the faces than the ball's distance: the centre of a ball resting on a face lies within its
distance of it).

Where it meets the other end ball first, or meets its end against the way the edge runs, the balls do not follow
along the path as along the edge, and the spine does not carry the rail either. Past the end balls the trace
goes on as far as the rail is carried past its ends, at the rate the spine's length between them runs to the
edge's range, and a step beyond, or until it leaves the domains.

On the class of rails below, the trace takes 0.2-10 ms a rail where the intersector's starting-point walk from
the middle ball took 3-23 ms and the projections of the middle ball onto the general offset surfaces up to 15 ms
more. The intersector's walk also split one branch into several lines (a cylinder and a torus: six), none of
which held all three balls; the trace follows the branch through them.

The station parameter is affine in the spine's arc length between the end balls, so the ball
at an end of the edge is the edge's own and the stations between are evenly spaced along the
centres' path. The march grows its first ball by the rail's own stations, which are the same
balls, keeps the existing side check on it, and then places it on the spine. The fit takes no
separate centre parameter on such a rail: its stations already have one.

### When the spine is used

Every open rail of that route is traced before it marches. A rail's own station at `t` is the ball whose plane
across its path passes through a point `P(t)` (the edge's point, or the centre across it). On a spine `c(s)`
with unit tangent `T(s)` by arc length `s`, that ball is where `h(s,t) = (c(s) - P(t)).T(s) = 0`, with `h_s = 1 +
(c - P).T'` and `h_t = -P'.T`, so `dt/ds = -h_s/h_t` (**derived**).

The station map is one to one exactly where `dt/ds > 0` along the zero set of `h` between the two end balls.
`h_s` vanishes where the focal line of the planes across the spine passes through `P` (the planes there turn
faster than `P` moves away from the spine), and `h_t` where `P` moves across the spine: at either sign change no
march by `t` passes, and the march stalls.

The zero set is traced from the first end ball by pseudo-arc-length continuation (Allgower and Georg, *Numerical
Continuation Methods: An Introduction*, Springer, 1990, chapter 2), with `t` scaled to length by the spine's
length per unit of the edge's range; a step is kept where the corrector converges near its prediction and its
tangent turns by less than 22.5 degrees, otherwise it is halved (a step policy, not an existence bound).

The rail walks the spine where the traced map turns back in either variable, and keeps its own stations where
the trace reaches the last end ball monotonically or cannot be read. There is no threshold: a sign. A closed
circuit keeps its seed through the edge's points and its affine parameter correspondence unchanged.

The spine itself is the offset intersection of Choi and Ju, *Constant-radius blending in surface
modelling*, Computer-Aided Design 21(4), 1989, pp. 213-220, traced as in Patrikalakis, Maekawa
and Cho (above). Arc length is not a new acceptance argument: every ball is still corrected on
both original offsets and the fitted sheet faces the unchanged contact, crossing, attachment
and solid checks. In particular a spine that bends more tightly than the ball (`kappa r > 1`)
still meets the interior-crossing check, which refuses the self-crossing envelope: that
trimming is not built here.

On the groove of a torus (radii 9 and 4) fused with a B-spline surface skinned through three
ovals (the longest concave edge and its tangent chain, distance 1) the rail's own station map turns back
at s = 3.78 mm (t = 0.771) and the centre route's at t = 0.448. Along the spine the stations
are evenly spaced (13.7 mm between the end balls) and the spine's curvature peaks at 1.35/mm
near t = 0.52. The conic of rho 0.3 then delivers single and chain selections (13 faces, every
sampled contact within 3.3e-8 mm of distance 1 from both supports); the round and the G2
section meet the interior-crossing check there, since their sections reach where
`1 - kappa r cos a` changes sign.

### Radius-law extrema

Positivity and the scale used by the tolerances now come from outward-rounded enclosures of every cubic span,
rather than uniform samples. Writing a span over x in [0,1], with h its parameter width, endpoints y0,y1 and
stored second derivatives m0,m1, its Bernstein controls are derived here as b0=y0,
b1=y0+(y1-y0-h^2(2m0+m1)/6)/3, b2=y1-(y1-y0+h^2(m0+2m1)/6)/3, b3=y1. The derivative controls are 3(b1-b0),
3(b2-b1), 3(b3-b2). Their interval signs certify monotone subspans, whose extrema are endpoints.

Otherwise de Casteljau subdivision isolates every possible derivative root, including multiple roots, until the
coefficient intervals overlap or the parameter interval has no representable midpoint. The remaining coefficient
hull encloses its stationary values. Thus every span endpoint and stationary region contributes to both bounds;
no fixed sampling width can miss an excursion. Nonfinite bounds or an unresolved positive sign report a
radius-law failure before any ball construction.

The range enclosure follows the nonnegative partition of unity of the Bernstein basis; the derivative and
subdivision argument above is derived here. Read: Fuenfzig, Michelucci and Foufou, [*Polytope-Based Computation
of Polynomial
Ranges*](https://cfuenfzig.lima-city.de/Publications/PolytopeBasedComputationOfPolynomialRanges_CAGD2012.pdf)
(2012; preprint 2011), sections 3.1-3.2: interval extensions and Bernstein convex-hull bounds. The existing
interval primitives enclose arithmetic on stored binary64 inputs under round-to-nearest without fast-math.

These bounds concern the stored spline, not uncertainty in the input data or the linear solve. Native checks
cover a narrow negative natural-spline span, 18 quadratic extrema across independent span widths and signs,
three multiple derivative roots and a constant law.

### Steep and flat wedges

In a wedge of opening psi between two planes the ball's centre lies r / sin(psi / 2) from the edge (the centre
is r from both planes, on the bisector), and its contacts r cot(psi / 2) from it. Past an edge's end the blend
is carried on while the ball goes on with the rail; the carry stopped where the centre was more than 4 r from
the rail, which a ball in a wedge sharper than 2 asin(1/4) = 29 degrees is already at the end, so its contacts
stopped short of the faces that stop the edge.

The bound is now that many times the ball's own distance at the end over sqrt 2 where the faces meet more
sharply than at a right angle (unchanged at a right angle and wider), so it no longer depends on the opening.

Measured over straight edges between two planes, a plane and a cylinder, and two cylinders, convex and concave:

| Measure | Result |
| --- | --- |
| Blends | 5040: openings 5 to 175 degrees in steps of 5; at rest, turned and scaled by 1e-3 and 1e3; Round, Conic 0.7 and Chamfer |
| Families | the plane-plane family was used to find the cause; the cylinder families were held out |
| Refused at openings below 30 degrees | 516 before, 12 now |
| The 12 still refused | chamfers of setback 1 at a concave opening of 170 or 175 degrees beside a cylinder of radius 10, whose chord between the feet runs through the cylinder (refused by name, as before) |
| Volume removed or added by the planes | r^2 (cot(psi / 2) - (pi - psi) / 2) per unit length (Round) and d^2 sin(psi) / 2 (Chamfer), to 3e-12 of the body |

Where the opening turns along the edge (the rim of a hole of radius R leaning a from upright through a plate)
the centres lie on an ellipse of semi-axes R + r and (R + r) / cos a, whose radius of curvature at the acute end
is (R + r) cos a.

The rolling-ball surface (the envelope of the balls, a canal surface) is regular where 1 - kappa r cos(alpha) is
non-zero, alpha the angle of a point of the section from the path's principal normal (derived here: along the
path the surface c + r (cos alpha N + sin alpha B) moves by |c'| (1 - kappa r cos alpha) along the path's
tangent, its other terms lie along the section; it is the regularity condition of a pipe surface in Maekawa,
Patrikalakis, Sakkalis and Yu, CAGD 15, 1998, not read here); on the rim the arc between the contacts passes
through that normal at the acute end, so the surface crosses itself once r >= (R + r) cos a: from a = 70.5
degrees for r = R / 2, 78.5 for r = R / 4.

The contacts themselves do not fold there (their speed along the path is s'^2 + R (R + r) |u'|^2 over |c'|,
positive). Just short of the crossing the fit missed its tolerance (lean 75 for r = R / 4): there the section
turns by up to 165 degrees, and a single rational quadratic's middle weight cos(turn / 2) falls to 0.13 with its
apex 7.6 r off the arc, which the fit cannot follow as it moves.

A conic section that turns by more than 120 degrees anywhere along the edge (the faces' normals at its points;
along a whole chain) is therefore made of two spans, split at its middle by de Casteljau on the homogeneous
poles, the same curve (the 120 degrees, where the weight is a half, is a judgement; the 5040 straight wedges are
unchanged to 1e-9 of their volume). Where the surface does cross itself, the fitted sheet's sections run back
against the way the balls roll (S_u . (n0 x n1) < 0, sampled 65 x 7 places): refused by name.

Measured (the same 16 rims): the blend makes all four rims at 70 degrees and the two with r / (R + r) < cos 75
at 75 as the balls' envelope (every point r from the centres' ellipse to 1e-6), and refuses the crossing ones,
exactly the criterion above; trimming them is open.

### A size that changes along the edge

(`Request::radii`, `blendRange`). The size r(t) is the cubic spline through the given places with no curvature
at its ends (the natural spline: the peak 1, 3, 1 on a box edge is r = 1 + 6 s - 8 s^3 on its first half),
periodic round an edge that closes; and every section is the section of one size, r(t), in the ball's own plane
across the edge, the ball of that size resting on both faces with its feet there (the sections are circles there
to 1e-12).

The blend is G1 to both faces: at a foot its tangent plane holds the section's tangent, normal to the face's
normal in that plane, and the contact's tangent, which lies in the face. Past the edge's ends the spline goes on
along its end tangents (its second derivative is zero there, so it stays twice differentiable), and the ball is
carried there only while its size stays within half and twice its size at the end (a judgement: a steep law
reaches nothing within twice the blend's size); the carry is twice the size at that end.

The places are breaks of the fit (the spline is C2 there and the stations follow its value). Not taken: the
rolling-ball envelope of balls of radius r(t), whose characteristic circle lies in the plane (x - c).c' = -r r'
(the derivative of |x - c|^2 = r^2; the feet lie on it since n.c' = r'); it was built first and differs from the
chosen section by 0.06 mm^3 on a box edge with a linear law, and with the cubic of `Law_Interpol` (OpenCascade's
variable fillet's interpolation) by up to 2.4 mm^3 on a peak.

A chamfer's setbacks follow the same law; at an angle off a reference face the other setback is d sin a /
sin(theta + a), theta the opening between the faces' tangent planes at the edge's point (the triangle of the
edge's point and the feet). Distances off the two faces of a round, conic or G2 blend must be equal where the
size changes (refused by name otherwise); such edges add no tangent edges and take the plain construction (no
island, re-blend or mixed-convexity routes).

A range keeps the change the blend of the edge alone makes between the planes square to the edge at its limits,
as `filletRange` does for the round. Where the premise below holds, the limited change is built directly
(`directPremise`, `directRange`).

The premise: the blend is this module's own, of the one edge, whose two faces each change into one face; the
blend makes one new face (the sheet) that meets each of them in one edge (its contact); both limits lie strictly
inside the edge; every other edge the blend changed or made lies strictly outside the slab between the limits,
and each contact enters and leaves the slab once, through the two limits (sampled, as the range's other checks
are; not an interval certificate).

Within the slab the change is then the strip of each face between the edge and its contact and the sheet between
the contacts, so each limit plane meets it in three curves: across each face from the edge to the contact, and
across the sheet between the contacts.

The result is the body in which, in each face's wire, the edge's piece within the slab gives way to the curve
across the face, the contact between the limits (the blend's curve and both of its curves on the surfaces,
restricted in range only) and the curve back; the sheet between the limits enters with the side the blend gave
it; and each limit is closed by a planar cap bounded by its three curves, the cap using the sheet's section in
the opposite sense.

A plane meets a plane, cylinder, cone, sphere or torus in lines and conics, taken exact (OpenCascade's analytic
intersection) and their curves on the face and in the cap by projection, which is exact for them.

On other faces each sample of the section is one Newton step sequence across, in the face's parameters, onto the
plane; where the solved points lie on an iso curve of the face to rounding, that iso curve is the section;
otherwise the curve on the face and the curve in the plane are interpolated at the same parameters, refined
while they come closer, and the edge takes their measured distance as its tolerance.

Each new vertex takes the largest distance from it to the ends of its curves and of their curves on the
surfaces, at least 1e-7, and the construction is refused where a gap exceeds the tolerances its curves carry (a
construction that does not close). Every other face, edge and vertex is the body's own, and the faces go back
into the body's shell. Outside the premise (the stock builder, a limit at an end of the edge, a blend of several
new faces), the change is isolated by Boolean operations as follows.

Where the blend reports which faces it kept unchanged, those pairs cancel and only the remaining skins are
intersected, once, in one non-destructive intersection structure shared by both oriented differences (the
separation of intersection and building, Open CASCADE SAS, *Boolean Operations*, sections *Intersection Part*
and *Safe processing mode*).

The retained pieces of the source faces are substituted in place and the new sheets enter reversed; the pieces
are grouped by their shared topological edges, never by a tolerance merge or sewing, so the source surfaces and
the edges' curve representations are kept as they were. Each group is then a closed shell of the removed
(positive signed volume) or the added (negative) material. **Derived here:** that grouping gives whole outer
shells when both bodies are one solid bounded by one shell and share an unchanged face.

A cavity shell of O\N with a face on the old boundary would put points outside O inside a region enclosed by
O\N, so O would have a void; one with all faces on the new boundary would hold all of the connected body N, the
shared face included, which then lies on the old boundary inside that region: the same contradiction. N\O
follows with the roles exchanged. Bodies outside that domain, and the stock builder, intersect the whole bodies
instead. The stopping-plane clips and the final cut and fuse are Boolean operations.

Measured, direct against Boolean isolation on the generated family (box edge and cylinder rim, convex and in a
cavity, conic, G2 and chamfer, scales 0.01, 1 and 100, limits 0.17 and 0.63): the same face kinds, faces within 5.2e-8
both ways and volumes within 4e-12 of each other, except where isolation was wrong: a large G2 box in a
cavity it returned unblended (the direct change is the range's share of the whole blend's change to 3e-11, as at
the smaller scale, the sections being the same all along) and two G2 rims in a cavity whose trim curves ended 1.2e-6
past their vertices. Vertex tolerances fall to the real gaps: 1e-7 at scale 1 where isolation left up to 1.1e-4.

Measured (`BRepBlendExt_Test`):

| Case | Result |
| --- | --- |
| Rounds on a box edge (linear, steep, peak, valley), at rest, turned and at scale 1e-3 and 1e3 | on their circles to 1e-6; remove (1 - pi / 4) times the integral of r^2 to 1e-9 of the body |
| A cylinder's rim (periodic) and, held out, the ring of a sphere on a plane | on their circles to 1e-6 |
| Chamfers of a growing setback, also at 30 degrees off a reference face and within a range | remove the integral of their sections to 1e-9 |
| G2 with a variable size or a range | seams G1 to 0.01 degrees and G2 to 1e-4 / r |

### Why spatial stations stay

`filletStations` (`BRepFilletAPI_MakeFillet` with a `Law_Interpol` law through sizes at points in space, cut at
limit planes) stays beside the section law, which cannot express what it takes: sizes by length along a chain of
tangent edges or round a closed contour (the section law holds one law per edge), and limit planes that are not
square to a straight edge (its range ends at sections across the edge).

Where both apply, on a box edge, measured both ways from tessellation points to the other body's mesh
(deflection 2e-3): a constant size between two limit planes gives the same body (4.8e-12 mm, the same volume) at 0.55
to 0.8 of the CPU time; one limit plane 7.8e-5 mm, the same volume; sizes 1, 2, 3 at three stations 0.046
mm and 0.53 mm^3 apart (`Law_Interpol` is not the natural spline); sizes at four stations within two limits 4.5e-3
mm one way; a constant size round an oval rim of four Bezier edges 1.6e-3 mm, at 1.5 times the CPU time.
G1 and G2 of these pairs were not measured.

### Connected limited contours

Plane intersections of each contour edge partition its parameter range into intervals of constant halfspace
membership. The existing length tolerance shifts the classification plane by that distance; topology is judged
at the same tolerance as before. Intersections use OCCT's
[GeomAPI_IntCS](https://dev.opencascade.org/doc/refman/html/class_geom_a_p_i___int_c_s.html) (Open CASCADE 8.0.1,
2026; `Parameters` and `Segment`), which reports the curve parameters of intersection points. This
remains a numerical intersection solve, not an interval certificate.

Count retained runs explicitly; merge the first and last only on a closed contour. This connectivity argument is
derived here. Two transitions alone do not imply one run on an open contour: its ends can both lie inside.
Native negative controls use the parabolic rim y=10(2t-1)^2, at three stopping heights 1, 5 and 9; each leaves
two separated end intervals.

### Over a side narrower than the blend, and a pinched body

Where the contact on a face lies past the face's far edge along the whole edge, the face is used up: a round's
ball rests on the far edge instead (a cliff), any other section is cut by the face beyond the far edge (a
notch). A contact that lies on the far edge itself (a side exactly as high as the ball is wide) is taken the
same way, as the ball there rests on the face and its edge at once.

Round a closed edge (a cylinder's rim over its side) the notch's edge on the face beyond is a loop of its own,
which takes the far edge's place, and the blend's seam section ends on it.

Measured: the top rim of a cylinder of radius 10 whose side is 0.5 to 2 high, rounded or chamfered by 1 to 3 at
least as high, keeps three faces and removes pi times the integral of R^2 - x(z)^2 (the arc by its angle, the
bevel in closed form) to 1e-9 of the body (round) and to the length tolerance over the bevel (its edge on the
bottom is interpolated); the same rims (side 1 to 8, radius 1 and 3) agree with independent meshes within 6.5e-7
mm for every round and chamfer, and G2 at a side as high as the ball within 1e-3 mm^3 of its volume.

Where the side is lower than the ball the blend makes G2 too (the G2 section cut by the bottom face, its seams
G2 to 1e-14 / r). A hole tangent to a face beside the edge pinches the body: the hole meets the face in a line
along which four faces meet. No blend passes it and the blend is refused by name; with a gap of 0.25 to 2 every
round, chamfer and G2 lies within 1.3e-6 mm (G2 1e-4) of an independent mesh.

### A run back of no length

A dip in a contact's speed the way the balls roll is a run back only
where the stretch of contact between the zeros of that speed either side of it (bisected) is longer
than `together` (1e-6 of the size), its length summed over 16 chords. Where the fit's own parameter
stalls (a stationary point of the sheet's parameterisation) both contacts stop together and go on the
same way, and the sampled least is rounding on no length: on a torus beside a B-spline face, |dS/du|
at the middle section was 2.6e-2, 2.6e-3 and 4.5e-8 at du = 2e-3, 2e-4 and 0, both contact speeds
positive on either side, run back 8.5e-14 mm. Those requests now refuse for their own reasons further
on. The 16 chords are a judgement (a chord sum is below the arc's length).

### Where a contact folds on any face

A contact folds where the blend's size passes the radius of the face's curve it follows (a setback r on a
cylinder's wall beside a rim leaning a degrees off the axis folds once r passes the rim's geodesic radius of
curvature on the unrolled wall, about R / tan a): its run back lies between two cusps, and the contact crosses
itself at X. The crossing is found in the face's parameters whatever its surface (formerly planes only), and the
blend is cut at the two sections through X. A round's loop is patched by its ball at X. A chamfer is the outer
envelope of its sections (on the cap, the contact trimmed at X: the distance set), so its two sheets meet in
their crease, S(ra, va) = S(rb, vb) on the one sheet, traced from X as the two sections close to the sheet's
swallowtail point V; from V to the other contact the two faces meet along one section, smoothly.

No patch. Where the loop's tail leaves the face over the blended edge itself, the face-crossing scan goes on
past the loop, which is cut away. A setback's side of the edge is read from both faces' normals at the edge's
point, not at its feet (a foot a setback round a curved face has its normal turned by the setback's angle, which
where the faces meet nearly flat or nearly closed turns the feet's direction over), and a chamfer contact's
speed is read along the edge, the normal of its section plane.

Measured: the slanted rims of a 4 or 8 thick plate, leaning 70 to 85 degrees, hole 1 or 2, setback 0.3 to 1, are
chamfered in all 16 cases: every contact point lies at distance r from the rim within 8.6e-8 mm and no blend
point lies beyond r; the removed volume matches the union of the section triangles (Monte Carlo, 2e5 points)
within 2.5 standard errors (0.6 to 0.8%). The planar-cap chamfers of the trimmed-outline corpus keep their
distance-field volumes.

### Where a round crosses itself between its contacts

The pipe c + r (cos a N + sin a B) about the
balls' centre path moves at |c'| (1 - kappa r cos a), so where the path bends more tightly than the
ball (kappa r > 1) the sheet folds over the middle of its sections while both contacts roll on: the
cuspidal edges 1 - kappa r cos a = 0 enclose the fold. A round is the boundary of the balls' union,
the distance set at r from the centre path, as a chamfer is the outer envelope of its sections: the
plane offset trimmed at its self-intersection is the analogue (Farouki and Neff, *Analytic
properties of plane offset curves*, CAGD 7 (1990), sections 1 and 4). Over the fold the sheet before
it meets the sheet after it in a crease S(ra, va) = S(rb, vb), ra < rb, and every section between
lies inside the union. The crease's two places run round one loop in the sheet's parameters, closed
at its two ends, the sheet's swallowtail points, where the sections before and after are one and
which lie inside the sheet; from each end to the contact beside it the two sheets are one section,
so the blend is cut along it. The edge's piece becomes two (a fold whose joint carries the crease),
each on its own copy of the surface, and each face's end is three edges: the section through one
end, the crease (G0), the section through the other. Further cuts along smooth sections are not
needed: the two faces are the fewest the crease allows.

Construction (Blend::sheets through checkCrossing, crossings; Attach::crossingFold). The fold is
read where S_u turns against the rolling direction on 65 x 7 samples. On the section line v = v0 of
its deepest sample, the sheet before the fold meets the sheet after it (GeomAPI_IntCS) in reaches
about the fold that double from its width, as a cap's loop is found; the narrowest meeting seeds the
crease. From there it is traced towards the end where the gap rb - ra closes, the gap held and cut
by a share each step (the chamfer's crease, Attach::creaseFold), to where the two sections lie
within the length tolerance of each other; towards the other end by pseudo-arclength steps (the
null vector of the 3 x 4 Jacobian by its signed minors; Allgower and Georg, *Introduction to
Numerical Continuation Methods*, 2003, chapters 2 and 6) past the widest gap, then closed the same
way. Near an end the held-gap system's least singular values fall with the gap; they are kept down
to rounding (the SVD's default cut-off at 1e-6 of the largest left those directions to the step's
prediction: the crease's end moved with the step size by 3e-5 of the sheet's u, and a face took in
4e-4 mm of the fold). The crease's three curves share one parameter, the length along its places
before the fold in the sheet's parameters (scaled by the sheet's speeds there): the loop is smooth
through its ends, unlike the crease's length in space, in which the places go as its square root.
Where a place crosses a knot line at which the sheet is C2 or less (the fit's C1 knots at support
breaks, the C0 join of a section's two arcs) the crease is solved there with that parameter held,
and the curves are interpolated piecewise between such places and joined: a cubic is
fourth-order accurate only on smooth data (de Boor, *A Practical Guide to Splines*, 2001, chapter
V). A span whose curves part from the sheet by more than the length tolerance is halved by a place
solved on the plane across its chord, as interpolated() halves a parameter curve's spans. Refused,
with the earlier message, where the crease leaves the piece, a place crosses the section through an
end (the faces cannot be cut there), or the curves would need more than 4096 places. Other laws keep
the refusal: their sections are not the balls', so the union does not define them, and a crease
inside a G2 blend would break the curvature continuity it is asked for.

Unproven choices: a place found between two must lie within half their chord of its middle and keep
a gap of at least a quarter of the wider of theirs; the first arclength step is an eighth of the
seed's gap, growing by 1.5 up to it; a step turning the direction by more than 45 degrees is halved.
Measured on a torus fused with a B-spline surface skinned through three ovals, along their concave groove (kappa r up to 1.35 at
distance 1): rounds at distance 0.98, 1 and 1.02, single edge and tangent chain, are made with 14
faces where they were refused; every tessellation point of their blend faces (4951 to 5648) lies at
distance r from the centre path within 5.0e-7 mm (the tessellation's own error) and none nearer, and
every contact's ball is r from both supports within 2.8e-8 mm. In a box of side 0.5 about the
crease, clear of the supports, the result's volume matches the distance set's, integrated by
columns: they differ by 1.8e-4 and 9.2e-5 of it at 60 x 60 and 100 x 100 columns, which extrapolates
to agreement. The rounds at 0.9 and 0.95 are unchanged.

### Further constructions

Three constructions go beyond the common cases: a chamfer whose
distances vary along the edge, G2 where a narrow side is lower than the ball (h < r), and a chamfer larger than a
round it runs over. Each is held to the same definitions and checks as the others.

### Unproven choices

Every marching limit (step counts, stall windows, Newton tolerances), the carry's stray bound
(4 r, or 4 / sqrt 2 times the ball's distance at the end in a sharper wedge), the radius/curvature
refusal margins, the choice between arcs and spline sections at corners, and the fit's pole and
stagnation budgets. There is no published bound for a G2 conic or quintic section as a function of the support's
curvature [searched; none found]; curvature-continuous requests therefore refuse in this tree until a caller
supplies a policy. Where a result is delivered it is measured again from outside (seam gaps, tangent angle,
curvature difference) and refused if it misses; that measurement is sampled, not a global certificate.
The separation certificate's split rule (the larger pole box, across its longer control polygon) affects
only its cost.

### Stock counterpart

`BRepFilletAPI_MakeFillet` (ChFi3d). See the table below.

Corner contact stations likewise compare the integral of rail speed between
their normalized parameters, rather than a raw parameter difference. The
arclength bounds the chord displacement and has physical length units; a
parameter interval need not. The existing agreement bound is unchanged.

### Local refinement without an iteration budget

A request's default `FitWork::Deferred` first runs the complete existing blend
pipeline with its historical sixteen-round fitting policy (`FitWork::Bounded`).
Only after all its constructions refuse geometrically does it repeat with
`FitWork::Extended`, where short converging intervals may refine further.
Requests made while reconstructing another blend inherit this policy. A host
with additional public constructions must request Bounded first, try those other
constructions, and only then request Deferred. Thus deferred work cannot replace
a successful existing construction with a different representation. If both
passes refuse, the original refusal is retained; cancellation propagates.
This selection rule is derived from ordered choice, independent of geometry.

The later passes (extended fits; refined contacts; then offset spines) each change one input of the request, and
each runs only where the pass before it read that input (`StageReads`). A bounded fit reads its policy where it
stops at its sixteenth round or goes on past a round beyond the eighth that inserts no knot: everywhere else the
extended fit takes the same steps. Attachment reads the refined contacts where they take another branch
(`Attach::refines`) or give a trimmed contact an owner.

Each place where offset spines would act records whether its spine is adopted (traced and its stations fold) or
its terminal contact could be continued; these are asked only once the pass has refused. Before the first read
of a changed input two passes compute the same thing from the same data, so a pass whose changed input was not
read repeats the one before it step by step and refuses where it did (derived, by induction along the
computation). Only the first refusal is reported, so skipping such a pass changes no result.

That refusal gives way to a later pass's only where it names a construction that failed (a station not found, a
fit not met) and the later pass measured the obstruction on its sheet (it crosses itself, a contact folds): on a
hole's rim through a plate leaning 80 to 85 degrees, the G2 stations fail where the balls' centre path bends at
kappa r = 1.15 to 3.83, which the offset spine's sheet measures as its crossing. Measured on 27 refusal rows of
torus, cylinder, extrusion and B-spline supports, every pass so predicted refused as before.

The existing interpolation, tolerance, continuity events, positive-weight check
and downstream result checks remain authoritative. The sixteen-round pass is
retained for construction selection, not as an existence or error bound.

Carl de Boor, *Spline Basics* (2001 manuscript for the *CAGD Handbook*),
[sections 14 and 17, theorem 17.2 and equation 17.3](https://ftp.cs.wisc.edu/Approx/cagdhand.pdf),
gives knot insertion and the Schoenberg-Whitney existence conditions for point
interpolation. This supports collocation under those conditions, not convergence
of adaptive rational or terminal Hermite fitting or a whole-domain certificate
from sampled data.

Termination argument (derived): scheduled break discovery occurs at round 2,
regime discovery at round 3, and bend-only checking expires before round 8 (or
sooner at its existing bend-round budget). In the extended pass an insertion
stall advances through these finite scheduled transitions. After round 8,
a continued ordinary refinement must add at least one distinct representable
knot. Otherwise the fit stops. Positive knot multiplicities imply that increasing
distinct knot count increases coefficient count; the unchanged 4,000-coefficient
budget therefore bounds continuing rounds. Event rounds may replace nearby
optional knots, but occur only twice. Cancellation and stagnation still apply.
This proves finite work, not convergence or existence for arbitrary contact data.
No result bound is relaxed.

### Located continuity events in a fitted blend

A newly located break is inserted at its own parameter. Nearby optional
simple fit knots may yield to it; a previous continuity event is retained.
The existing knot-neighbourhood size only selects optional fit knots to remove,
not a displacement allowance for an event. **Derived:** moving a true break
by delta changes the physical station by the integral of rail speed over
delta, so a normalized parameter threshold alone cannot bound that movement.
The final fitted-surface and contact checks remain unchanged. Very close
but distinct events can still require an ill-conditioned or costly fit;
this change does not identify them as the same event.

### Construction boundary roles

Construction roles are now retained separately in `Blended::joins`: contacts,
used-support joins, joints, support breaks, folds, rolls, corners, sharp mitres,
stopping ends, trimmed sides, ridges and split boundaries. Each carries its required
continuity order, independently of `BRep_Builder::Continuity`. Recording an order is
intent, not evidence that the delivered geometry attains it. Fold and roll contacts
require tangency, while a support break does not claim curvature continuity.

The inventory is filtered by final edge identity at attachment, and carried through
subsequent modifications in staged blends and Boolean restoration. This follows
OCCT, *Reference Manual* 8.0.1 (2026),
[`BRepTools_History`, detailed description and `Merge`](https://occt3d.com/dev/doc/refman/html/class_b_rep_tools___history.html):
modification, generation and removal are distinct relations. The edge-intent
propagation rule is derived here: use retained edges and modified edge descendants,
not unrelated generated shapes. The alternate face-envelope constructor has not yet
supplied this inventory; absence must not be interpreted as a smoothness claim.

For the geometric interpretation, Patrikalakis, Maekawa and Cho (2009),
[section 8.1.3.5, Linkage Curve Theorem](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node156.html),
assumes tangent-plane continuity along a first-order continuous boundary before
curvature agreement in a transverse direction establishes second-order smoothness.
It does not classify construction intent from the measured output angle. Native
role checks cover twelve combinations of three scales, two noncircular sections,
and single stopped edges versus closed mitred loops.

### Compatibility before a second-order join

Jorg Peters, *Geometric Continuity* (2001), section 2, defines continuity by regular reparameterization, rather
than equality of arbitrary parameter speeds: <https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf>.
**Derived here by the chain rule:** orient a common seam as `A(0,s)=B(0,s)`, and write the first transverse
reparameterization derivatives as `alpha,beta` and the second as `gamma,delta`. Then `B_t=alpha A_u+beta A_v`
and `B_tt=alpha^2 A_uu+2 alpha beta A_uv+beta^2 A_vv+gamma A_u+delta A_v`.

The mixed derivative follows by differentiating the first equation along s. The parameter map must remain
nonsingular; in these coordinates alpha is nonzero, with its sign fixed by the chosen orientation. Contact and
seam data meeting at a corner must be derivatives of the same map. A quintic spline with interior multiplicity
at most three is C2, but this alone imposes neither those support constraints nor regularity or fidelity to a
prescribed section family.

Salvi, Varady and Rockwood, *Transfinite surface interpolation over specific
curvenet configurations* (2013), section 3.2, equation (3.1), uses parabolic
ribbons `P(s)+d T(s)+d^2 C(s)/2` for curvature continuity; sections 3.2 and 4
distinguish compatible corner data from rational correction terms:
<https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf>.
This supplies a construction framework, not a convergence or error theorem for
the current section interpolation. In particular a fit's positional error does
not bound its second derivatives.

### Derived negative control

the graphs `(u,v,0)` and
`(u,v,.032 u^2 v^3(1-v)^3)` share a tangent plane along u=0, but their normal
curvatures across that seam differ by .001 at v=.5. Their contact ribbons at
v=0 and v=1 agree through order two. Reparameterization cannot remove this
geometric curvature jump while preserving the entire piecewise surface.
A native tensor-product polynomial evaluation gives normal jump zero,
curvature jump .001 and maximum sampled contact-jet roundoff 3.13e-17.
The compatible planar target gives zero curvature jump. This control validates
the distinction; it is not a classification of arbitrary blend requests.

### A remaining degree of freedom in symmetric quintic sections

For contacts A and B and tangent intersection T, let d=1-rho. Consider the
polynomial quintic basis with knots (0,.5,1), multiplicities (6,2,6), endpoint
poles A,B, adjacent poles `T+(1-.3 rho)(A-T)` and its symmetric counterpart,
and second poles `T+q(A-T)` and `T+q(B-T)`. Write the central pair as
`T+z(A-T)+b(B-T)` and its symmetric counterpart. **Derived here:** midpoint
interpolation gives `q+3(z+b)=4d`; the chosen reflected-pole relation
`b=(q-2d)/120` then fixes z, but leaves q free. Endpoint zero curvature,
symmetry and midpoint position therefore do not uniquely specify this family.
The reflected-pole relation is a construction choice, not a continuity theorem.

The endpoint derivative identities are `C'(0)=10(P1-P0)` and
`C''(0)=80(P2-2P1+P0)`. They follow from the standard B-spline derivative
recurrence: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer
Aided Design and Manufacturing* (2009), section 1.4.2, equations (1.65)-(1.68):
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html>.
Collinearity of the first three poles gives zero endpoint curvature; a curved
support requires different endpoint jet data.

Lian Fang, *G3 approximation of conic sections by quintic polynomial curves*
(1999), section 2, equations (1)-(2), constructs a quintic Hermite approximant
to a rational quadratic with endpoint G3 conditions:
<https://numerika.fmf.uni-lj.si/bib2html/e-documents/clanki/FangG3approximationofconicsectionsbyquinticpolynomialcurves.pdf>.
Those endpoint conditions differ from the zero-curvature contacts above, so
that published construction does not determine the free q here. Equating
midpoint curvature to a rational quadratic with the same contacts, tangent
intersection and midpoint is a separate, derived condition; it can have an
inadmissible root outside the contact triangle. No such rule is adopted.
A general generator and its existence conditions remain open.

### Internal curvature joins and boundary-ribbon compatibility

A quintic spline knot of multiplicity four permits a second-derivative jump;
retaining such a knot from an approximate intersection curve does not establish
that the underlying supports require a curvature break. Conversely, a curvature-
driven section law can lose rail smoothness when it crosses a merely C2 support
knot, although the support's curvature itself is continuous. Inspect the raw
section fit, joined sheet and retained trim separately before assigning a cause.
A positional knot-removal tolerance alone is not a curvature or contact bound.

For two regular polynomial patches with identical position and first derivatives
at a common seam, write D = S_u, T = S_v, N = D cross T and
J = S_uu(right) - S_uu(left). The mixed and along-seam second derivatives agree.
Equality of the second fundamental forms therefore requires N dot J = 0, not
J = 0. In a common orthonormal tangent frame the operator-norm difference is
abs(n dot J)/(E - F*F/G), where n = N/abs(N), E = D dot D,
F = D dot T and G = T dot T. This reduction is derived here from the fundamental
forms; it requires G > 0 and E*G - F*F > 0. It does not apply unchanged to a seam
whose first derivatives differ.

A possible fixed-space correction preserves position and first jets at every
seam and solves this scalar normal equation. With polynomial transverse degree
q, N has degree at most 2*q-1 and J degree q, so their scalar product has degree
at most 3*q-1 on each transverse span. Matching its Bernstein coefficients
provides an identity, whereas collocation only supplies sampled evidence.
Solvability still requires the right-hand side to belong to the constrained
coefficient map's range; more coefficients alone are not an existence proof.
Rational patches require the corresponding positive-denominator treatment.

Preserving both endpoint two-jets fixes the first and last three transverse
control rows of a clamped patch (after common knot refinement). A correction
vanishing to order two there cannot repair an existing nonzero endpoint
normal-curvature mismatch. Construct compatible, support-derived ideal ribbons
before freezing boundary constraints: freezing inaccurate fitted ribbons also
freezes their error. A merely G1 seam requires the full reparameterisation
conditions first. These are derived necessary conditions and a solver design,
not a new construction or an acceptance certificate.

Sources: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing*
(2009), section 3.3, equations (3.26)-(3.28), [normal
curvature](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node29.html), and section 8.1.3.5, [linkage
theorem](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node156.html); Peters, *Geometric Continuity*
(2001), sections 2-3, [author manuscript](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf);
Salvi, Varady and Rockwood, *Transfinite surface interpolation over specific curve-net configurations* (2013),
sections 3.2-3.3, equation (3.1), [parabolic ribbons and corner
compatibility](https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf). These sources establish the
geometric continuity and ribbon compatibility framework; the coefficient reduction and proposed correction above
are derived.

#### A common rail coordinate and support collars

An approximate intersection curve's stored continuity is not a regularity
requirement of its supports. A common rail coordinate can be inherited from the
carried curve, with affine changes on its pieces; it need not be arclength.
Lowering knot multiplicities in a positional section fit alone does not impose
support contact, repair a later first-order coupling, or prove convergence within
an iteration budget. A continuous target with a first-derivative break can need
linear, rather than high-order, local refinement. This last observation is a
limitation of applying smooth-target interpolation estimates, not a general
convergence theorem.

For a support parametrisation F and a contact parameter curve q, write a collar
as P(s) + d T(s) + d*d A(s)/2. Compatible support data have the form
P = F(q), T = DF(q) a, A = D2F(q)[a,a] + DF(q) b. The tangential fields a and b
are reparametrisation freedoms, not extra normal curvature. Differentiation
gives n dot P'' = II(P',P'), n dot T' = II(P',T), and
n dot A = II(T,T). These identities are derived here, with II evaluated on
physical tangent vectors and a regular support. The relevant traces of a C2
sheet have P in C2, T in C1 and A in C0. Fitting all coefficient functions in
the same derivative norm unnecessarily strengthens these trace conditions.
Neither these identities nor independent knot-wise constraints establish
contact between knots or compatibility at a shared corner.

A local diagnostic construction can share coefficient value, first derivative and second derivative at every
rail station. On a span of length h the quintic Bezier coefficients are P0, P0+h D0/5, P0+2 h D0/5+h*h A0/20,
P1-2 h D1/5+h*h A1/20, P1-h D1/5, P1. Shared jets imply C2 in exact arithmetic, even when redundant knots are
stored. Floating-point cancellation on small spans still needs a physical curvature check. For example, the
endpoint second derivative is 20*(B2-2*B1+B0)/(h*h).

If each stored coefficient component has absolute error at most e, its second-derivative error is at most
80*e/(h*h), before evaluator rounding; two independent sides can differ by 160*e/(h*h). These are derived
componentwise bounds, not curvature bounds: the surface metric and normal also enter curvature. A genuinely C2
knot representation avoids independent endpoint storage, but cannot supply missing support constraints.

In a common polynomial tensor Bezier basis, the maximum norm of corresponding control-point differences bounds
same-parameter displacement by the convex-hull property. It bounds neither normal differences nor second
fundamental forms. Knot insertion rounding must also be enclosed before calling this a numerical certificate.

After changing a sheet, build its contact curves and stopping rims from that
sheet. A projection fit's error in support parameters is not a world-space
error bound when the support metric is poorly scaled. Check the composed
support curve against the actual sheet in physical units, including endpoints;
sampled checks remain sampled evidence. Exact one-sided jets where rail knots
reach a retained contact can reveal errors missed by interior edge samples.
Preserving a positional proposal allowance cannot replace the original contact,
section, topology, regularity and internal-curvature acceptance tests.

Sources: Hartmann, *Geometry and Algorithms for Computer Aided Design* (October 2003), sections 11.3.2 and 11.4,
pages 133-135, [common offset-intersection
parametrisation](https://www2.mathematik.tu-darmstadt.de/~ehartmann/cdgen.pdf); Salvi, Varady and Rockwood
(2013), section 3.2, equations (3.1) and (3.3), [parabolic support
ribbons](https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf); Patrikalakis, Maekawa and Cho (2009),
section 1.4.3, [knot removal and
continuity](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node18.html).

These sources supply constructions under their stated regularity and compatibility conditions; they do not prove
existence of an arbitrary support-constrained approximation in a fixed tensor space. The trace identities,
Hermite diagnostic and preservation qualifications above are derived here.

#### Normal-acceleration corrections and numerical evidence

For a regular seam whose tangent planes agree, let t be its unit tangent and
c = n cross t. A correction delta*n to S_uu changes II(c,c) by
`delta/(S_u dot c)^2`, with the seam position and first derivatives held fixed.
This identity is derived from the inverse surface metric. Matching this scalar
alone is sufficient only if the along-seam and mixed second fundamental forms
already agree. First-order coupling of independently fitted pieces does not
establish those hypotheses; check the full form in one physical tangent frame.

A derived compact acceleration function on a one-sided interval of width h is `h*h*x*x*x*(1-x)*(1-x)/2`, where x
runs from zero at its outer boundary to one at the seam. It vanishes through order two at the outer boundary,
has zero value and first derivative at the seam, and has second derivative one there. Its maximum is
`54*h*h/3125`. Thus a coefficient field bounded by Q displaces the surface by at most `54*h*h*Q/3125` in exact
arithmetic. This is a positional proposal bound, not a curvature or support-contact certificate.

Transverse collocation of the field does not establish an identity between samples. Successive knot insertion
and world-coordinate coefficient rounding can destroy the intended second derivative on very short spans. Adding
more transverse sites cannot repair that independent longitudinal conditioning problem.

Use independent, well-conditioned derivatives before judging a construction. For a degree-p spline with expanded
knots U, its derivative coefficients are `p*(P[i+1]-P[i])/(U[i+p+1]-U[i+1])`, on the knot vector with the first
and last knots removed. Apply the formula in each tensor direction and to homogeneous coordinates for rational
surfaces, then apply the quotient rule. Translating positions by a fixed origin before evaluation leaves all
physical derivatives unchanged in exact arithmetic.

Form first homogeneous differences directly as `w_b*(P_b-P_a)+(P_a-O)*(w_b-w_a)` for a fixed origin O (derived
identity); subtracting the origin from each pole before differencing can itself discard the small difference on
a short span. Higher differences then act on derivative coefficients. Summing differentiated basis functions
against large world positions can instead cancel terms proportional to the inverse square of a short span.

Derivative control nets reduce evaluator cancellation; they do not recover information already lost in stored
coefficients. Increased-precision evaluation of the same stored net distinguishes these two sources of error.

A positional fit across an entire carried rail is likewise insufficient: making
all internal knots structurally C2 can approximate a tangential parameter kink
using very small spans while introducing severe support curvature. A structural
alternative is to solve compatible contact ribbons and their parameter freedoms
jointly, controlling the required derivative norms and the longitudinal
conditioning before fixing a tensor partition. This is a design direction,
not an existence claim for arbitrary boundary data.

Sources: Nowacki, Kaklis and Weber, *Curve mesh fairing and GC2 surface interpolation* (1992), sections 3.2-3.3,
equations (26)-(42), [quintic Hermite interpolation and
compatibility](https://www.numdam.org/item/M2AN_1992__26_1_113_0.pdf); Peters, *Geometric Continuity* (2001),
sections 2-3, [reparameterisation conditions](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf);
Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009), section 1.4.2,
equations (1.65)-(1.68), [derivative control
coefficients](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html), and section 1.4.4,
[tensor-product surfaces](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node19.html).

The compact correction, its maximum and the conditioning conclusions above are derived here; the cited sources
supply the interpolation and continuity framework, not a proof of convergence of this correction.

### Contact ownership on a joined blend sheet

A final construction, attempted only after the previous bounded and extended pipelines refuse, refines the
**final joined sheet's** ownership intervals. An earlier rail split is insufficient when fitting and joining
change the contact chart. Isolated original-contact/beyond-surface intersections are retained only when they
project onto the actual shared boundary edge within the existing tolerance. They partition increasing active
rail intervals; original-face IN/ON intervals retain that owner, other intervals retain the existing beyond-face
owner.

Adjacent pieces share a witnessed Split boundary. The 3D sheet is unchanged. Copied 2D fitted contact curves are
restricted by `Segment` at interior events, retaining the original outer continuation tails. A bounded notch
crossing one far edge twice retains the two outside edge intervals together rather than applying incompatible
sequential endpoint cuts.

The common-refinement and interval-set invariants are derived here. Carl de
Boor, *Spline Basics* (2001), section 14, establishes exact knot insertion;
restriction therefore does not approximate the contact curves
([chapter](https://ftp.cs.wisc.edu/Approx/cagdhand.pdf)). Sections 6, 10 and 17
also explain the nonnegative partition of unity, convex-hull bound and
interpolation existence conditions used when diagnosing missing endpoints.
Jaiswal, Ghosal and Gurumoorthy, *An Efficient Algorithm for Circular Blending
of Edges* (2001), introduction, page 1, distinguish blend-sheet generation
from trimming and topological integration; their local offset construction
does not supply this arrangement
([paper](https://mecheng.iisc.ac.in/~asitava/sharad.pdf)).

This is a numerical construction for finitely many isolated boundary events on a regular contact strip, not a
certified enumeration of every possible contact component. Overlap segments and coincident events on distinct
edges are refused, and decreasing active charts are not refined. Midpoint ownership classification presupposes
that the witnessed boundary events delimit the relevant owner changes; the unchanged final topology, tolerance,
continuity and collision checks remain mandatory.

Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009), section 5.8.2.3,
equations 5.100--5.105, explain why marching requires starting points for every component and why
tangencies require additional treatment
([chapter](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html)).

The OCCT `GeomAPI_IntCS` interface separately exposes isolated points and tangential segments; a zero point
count alone is not a complete empty-set oracle
([reference](https://occt3d.com/dev/doc/refman/html/class_geom_a_p_i___int_c_s.html)).

### Centre-path arc coordinates for a defective sheet

A rail-parameter knot can be a tangential acceleration break rather than a geometric curvature break. For a
regular centre curve c(t), put v=|c_t| and sigma(t)=integral v dt. The chain rule gives `p_sigma = p_t/v` and
`p_sigma_sigma = p_tt/v^2 - p_t (c_t.c_tt)/v^4` for each contact foot p(t).

This cancellation is **derived** from the arc-length relations in Patrikalakis, Maekawa and Cho, *Shape
Interrogation for Computer Aided Design and Manufacturing* (2009), section 2.1, equations (2.6)-(2.11), [arc
length and tangent vector](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node22.html). It removes only
tangential parameter acceleration; it cannot remove a geometric support break, a section-policy break, a
singular centre path or a varying radius law's break.

Peters, *Geometric Continuity* (2001), sections 2-3, [reparameterisation
conditions](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf), likewise distinguishes geometric
regularity from equality of parameter jets.

The experimental reconstruction keeps the existing result unless a sampled one-sided second-form comparison
exceeds the existing internal curvature bound at an inherited rail knot. For constant-size curvature-continuous
sections, centres are interpolated in the original edge coordinate, retaining its knots, and refined to the
station-solve length resolution. Numerical arc integration and inversion then define the fit coordinate; the
radius law and all station solves continue to use the original edge coordinate.

The experimental refit now keeps a checked reconstruction even if it needs more longitudinal coefficients.
Changing the fit coordinate does not preserve the complexity of its interpolant; coefficient count alone is not
an error estimate. Position and contact checks in `fit` remain mandatory. Additional coefficients are an
explicit representation cost, not evidence of a smaller geometric error or Pareto acceptance. The shared rail
carries both coordinate directions, so contact, section and attachment readers agree.

The sheet fit drops only edge-coordinate knots and retains support and policy events. Support-split rails retain
their previous construction: coupling their terminal jets requires a common group coordinate. Variable-size laws
and edge-supported contacts also retain their previous construction.

Position sampling of the centre carrier is not an interval bound on its length
or derivatives. The regularity argument is conditional on regular transverse
offset intersections and smooth support and section data; it does not establish
convergence of this sampled reconstruction. Final contact, internal curvature,
section-law, topology and preservation checks remain separate acceptance gates.

### Source-normal events of a centre carrier

The centre-based station plane uses both contact normals and normals evaluated
on the source edge's parameter curves. These are distinct paths through a
support's parameters. A contact can cross a support's finite spline boundary
before the source edge does. Retaining contact crossings alone therefore does
not retain every event of the centre function being interpolated.

For a regular station solve, the implicit-function theorem transports the regularity of its input data to the
centre. A support of class C^c has normals of class C^(c-1), so the centre carrier retains a quintic knot of
multiplicity `5-max(0,c-1)` at either path's crossing. A tangent continuation of a finite support is only C1 in
general and thus requires a C0 carrier knot. These regularity statements are **derived**, conditional on a
nonsingular station Jacobian.

Knot multiplicity and continuity follow the B-spline definition (`Geom_BSplineCurve`, Detailed Description). The
existing event locator is shared by the section and centre fits; the latter uses first-order support data rather
than the section law's second-order data. Connected pieces map their events into the delivering carrier's common
source coordinate.

During connected centre-carrier interpolation, station-plane origins outside a parameter curve's domain use the
existing adaptor's cached polynomial evaluation. The interpolator resolves positions to the station-solve
tolerance, substantially finer than the section fit. Ordinary station construction retains its original
parameter-curve evaluation. Each connected carrier starts from the delivering pieces' station seeds. Single-rail
centre fits retain their original event inventory and station cache.

Only a station interval extending beyond a source parameter curve's domain needs a separate cache evaluated with
the endpoint polynomial. The adaptor evaluates the same endpoint span, while avoiding repeated extrapolation
through signed de Boor coefficients. Horner evaluation's coefficient backward-error analysis is described by
Higham, *Accuracy and Stability of Numerical Algorithms*, second edition (2002), section 5.1, [Horner's
method](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch5).

It is not a universal forward-accuracy guarantee for extrapolation. Independent long-double de Boor comparisons
distinguish evaluation error from a geometric break. Neither the centre-fit length resolution nor any final
geometric bound is changed. Position sampling still does not certify the carrier's length or derivatives; final
contact, section, continuity and topology checks remain necessary.

### Continuity intent after limiting a blend

The optional join inventory of `blendRange` follows edges through each Boolean
modification used to isolate, clip and restore the material change. An edge
retains its contact order only along that edge's retained modifications. Faces
on the exact plane supplied to a clipping intersection create sharp closing
boundaries, recorded with stop role and order zero. Their curvature difference
from the blend is not a failed support contact. The geometry and the existing
solid acceptance checks are unchanged.

This propagation is derived from the Boolean history contract, not from measured
continuity or mutable edge flags. Open CASCADE Technology, *Boolean Operations*
(8.0, accessed 2026), section "History Information", distinguishes retained
modifications from new intersections generated by faces
([specification](https://www.occt3d.com/dev/doc/overview/html/specification__boolean_operations.html)).
The inventory does not certify continuity itself; consumers must still measure
the contacts and require coverage of the intended supports.

### Smooth section-policy selections

`BRepBlendExt_Smooth.hxx` provides compact C2 upper/lower rounding for policy
maxima and minima. On `|x| < 1` replace `|x|` by
`p(x)=3/8+3x^2/4-x^4/8`, retaining `|x|` outside. Values and the first two
derivatives match at both endpoints. Thus
`M(a,b,w)=(a+b+w*p((a-b)/w))/2` is C2 for positive width, is at least
`max(a,b)`, and exceeds it by at most `3w/16`; `-M(-a,-b,w)` rounds a minimum
downwards. These formulas and bounds are derived here by endpoint Hermite
matching and differentiating the excess. They apply to smooth scalar inputs;
a nonsmooth support-curvature field or a singular sweep is not repaired by them.

The G2 flat-support cutoff now activates a whole section profile with `H(t)=6t^5-15t^4+10t^3`, between one and
two fit-tolerance units of `abs(curvature)*radius^2`. Both endpoint derivatives through order two vanish.
Whole-profile interpolation is necessary when a host's flat and curved laws have different zero-curvature
limits. Contact points and tangent directions are shared by the two profiles; intermediate endpoint curvature is
eased too. The band width is a policy choice, not a new accuracy allowance.

Output contact, section, validity and collision checks remain independent requirements. Source: Ken Perlin, *Implementing
Improved Perlin Noise*, GPU Gems (2004), section 5.4, [quintic Hermite
interpolation](https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-5-implementing-improved-perlin-noise).
This source establishes the scalar interpolation regularity, not feasibility of arbitrary G2 blend data. A host
with C2 selections should return a constant section regime: changing labels would otherwise request unnecessary
C0 fit knots.

### Experimental terminal restriction continuation

The last recovery construction can continue a displaced support's terminal
contact along the identified two-offset branch. A local tangent carrier predicts
the next centre; the existing exact-offset Newton equations correct it. A sign
change against the stopping support is bisected and retained only when the exact
contact also projects onto its finite boundary arc within the unchanged bounds.
The endpoint records its support, arc, arc parameter and station parameter.
The unchanged attachment and final solid checks still decide delivery.

Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009), section 5.8.2.3,
equations (5.101)-(5.105), describe regular transverse-intersection marching and its component-seed
requirement ([chapter](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html)).

Isabelle Grignon, *ChFiDS_CommonPoint* (1993), `SetArc` and `SetParameter`, records separate restriction and
spine parameters in OCCT
([source](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_CommonPoint.hxx)).
The local predictor and signed-distance boundary search here are derived; the step size and finite exploration
count are numerical policies, not existence bounds. This WIP does not certify unsampled regularity, all boundary
roots, or consistent terminal face closure. The tested displaced contacts reach their restriction arcs; at this
stage complete solid construction still requires the boundary arrangement described below.

### Boundary intervals after two support crossings

In the final contact-arrangement construction, a boundary edge incident to a blended
chain can be crossed twice even when only one of its endpoints belongs to that chain.
The retained interval is bounded by the two crossing parameters. Trimming both times
from the chain endpoint loses one boundary vertex. Each cut now removes the endpoint
on its own side of the other crossing, using the ordered geometric parameters rather
than the parameter midpoint. This interval selection is **derived** for the isolated,
two-crossing arrangement; it is not a general arrangement algorithm for overlapping
contacts, tangencies or more crossings. Earlier constructions are unchanged.

The generated torus/oval-prism test checks single-edge and connected-chain requests
at two scales and two poses. It checks solid acceptance and, away from stopping
planes, the rolling-ball defect independently by distances to the input skin.
At larger radii the contacts can meet a stop plane at different longitudinal
stations, so replacing the trim by one iso-section would change the geometry.
The graph construction below handles regular planar branches across spline breaks.

Jacques Goussard, *BRepFeat_SplitShape* (1995), OCCT 6.8.0 edition (2014),
`Add(Compound, Face)`, requires geometrically connected splitting edges to share
their topological vertices
([source](https://github.com/Open-Cascade-SAS/OCCT/blob/V6_8_0/src/BRepFeat/BRepFeat_SplitShape.cdl)).
That boundary-incidence requirement motivates retaining both crossing vertices; it
supplies neither a complete topology proof nor a numerical error bound for this
construction. The existing final validity, closure and contact bounds still apply.

### Non-isoparametric planar stops across spline breaks

After the existing terminal intersection constructions refuse, the final contact
arrangement can trace a sheet/plane branch as a graph in the sheet's transverse
parameter. Scalar Newton correction solves the plane residual for the longitudinal
parameter. Both terminal witnesses must agree with the corrected branch within
the unchanged length tolerance. The branch is split at every transverse spline
knot and every intersection of a longitudinal knot line with the plane. Spatial
and parameter curves are then interpolated together on each smooth piece and
joined with their common parameter. The existing midpoint residual, solid validity,
closure and contact checks still decide acceptance. Earlier successful terminal
constructions are unchanged.

This graph chart and knot arrangement are **derived**. They require a regular branch that is a graph in the
chosen coordinate; Newton convergence, the finite refinement policy and midpoint checks are numerical tests, not
a global existence or approximation theorem. Tangencies, other components and branches leaving that chart remain
outside this construction.

Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
sections 5.8.1.1-5.8.1.2, equations (5.84)-(5.87), reduce a patch/plane intersection to a scalar parameter-space
curve and describe coordinate stepping with Newton correction, including the need to change charts near vertical
branches ([formulation](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node101.html),
[tracing](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node102.html)).

The generated two-radius, two-scale, two-pose single/chain test checks acceptance and the rolling-ball defect
against the input skin independently of the solver.

For an offset-intersection station route, centre-length fitting retains that
route's spline carrier, including its circuit seam and orientation. It does not
re-march the source edge. The carrier is segmented at circuit wraps and affinely
mapped to source coordinates before spanwise length integration. This coordinate
rule applies to every such G2 route, independently of a near-limit knot verdict.

### Carried-sheet corner boundaries

The affine three-sheet tensor is retained when all three solved contact crossings lie within the existing
positional allowance of the incident tangent planes spanned by the outgoing edge directions. This tests the
construction premise from geometry rather than the representation of the support. The coordinates, poles,
weights and boundary construction of that tensor are unchanged.

When a crossing is off its incident plane, the noncircular three-sheet construction uses the intersections of
the carried contacts as its vertices, and connects them by curves in each carried sheet's parameter domain.
These curves are approximated spatially and used as supported boundaries of the existing constrained fill. The
carried sheets are not changed. The two positional approximations share the existing allowance: if the endpoint
discrepancy is `g` and the total allowance is `B`, each receives `(B-g)/2`. The triangle inequality then
accounts for both approximations and the endpoint mismatch; this does not certify unsampled approximation
errors.

Tangent-continuous examples return valid closed bodies. The affine tensor
and the carried-sheet construction share one contact solve; the latter starts
from the same carried sheets and solved stations. Curvature-continuous
non-affine corners remain work in progress: the generic fill does not realize
the curvature-continuous ribbons within the existing bounds. Complete
derivative compatibility and a corner-corrected ribbon interpolant remain
unresolved.
No bound or interior acceptance gate has been relaxed.

Salvi and Varady, *Multi-sided Surfaces with Curvature Continuity* (2014),
sections 3.3 and 4.2--4.4, distinguish ribbon construction, compatible parameter
speeds, rational twist corrections and the blending order. Agreement of tangent
planes and second fundamental forms at the vertices does not imply agreement
of the higher mixed derivatives required by a regular tensor spline. The
current constrained-fill experiment is not an implementation of the paper's
corner-corrected Gregory sum, and supplies no such reproduction theorem.

### Stops reached past the edge's vertex

At a blend's end each contact is walked over its support from the edge's vertex
until it crosses an edge onto a face that is not a support. Two cases go past the
faces at the vertex.

#### A face of the support's own surface

A Boolean leaves a periodic surface's seam,
or a split, as an edge between two faces that are copies of one surface. Where the
contact crosses such an edge short of the face before it, it is on the same
surface, so the blend's sheet touches it there exactly as on the support. The walk
goes on over that face. The faces are taken as one where both give the same points
at the contact's parameters, from the crossing to the end of the carried contact
(16 samples to the length tolerance; this is a sampled test, not an identity proof).
The crossed edge is cut at the crossing. The contact is one edge on each face, joined
in one run along the blend's side, and each part keeps its own face's parameters,
shifted by whole periods where needed.

#### A face that stops only one contact

The blend ends on a face `F` only where the surface of `F` cuts the whole section. **Derived:** the end edge is
the meeting of the sheet with that surface, from a corner on one contact to a corner on the other. If the
surface meets only one contact where the sheet is carried, the meeting leaves the sheet before it reaches the
other contact, so `F` cannot end the blend. A stop surface carried past the vertex to one contact is therefore
taken only where it also meets the other contact.

A contact that crosses none of its face's edges where it is carried, and on whose face no other piece of the
blended edges lies beside the other support, ends on that face. (Where such a piece lies on the face, the ball
rolls on over the faces between to it, as before.) The face the other contact reaches then stops the blend, and
its surface is carried to the open contact. The open corner is where that contact first meets the surface, past
where it entered its face, and it must lie on that face.

The stop face grows to the blend: its new edges are its surface's meetings with the faces from its own vertex on
the first edge between them to the open corner. Each further edge between them is cut where the surface crosses
it. What lies under the blend goes: the swept edges, and the parts cut from those faces.

Measured on a plane meeting a crossed surface of revolution, concave, radii 0.2, 0.5
and 1, with circular, curvature-continuous and conic sections: all nine return valid
closed solids.
- Each round's blend face lies on its rolling balls to 1.3e-8 mm at 401 points per face.
  Each ball centre's distances to the plane and the wall are computed analytically.
- The added volume matches the integral of the sections' areas along the centre curve
  (tube formula, `1 - kappa N.x` weight) to 0.3 %, and a Monte Carlo union of sections
  within 0.7 standard errors.
- The corners at both ends agree with an independent display mesh to 1e-3 mm.
Before this, all nine were refused.
No bound was changed, and every earlier route keeps priority where it applies.

### Split corners of more than three conic blends

Where more than three conic blends meet at a vertex and every face there is a plane, each blend's
sheet is its section's extrusion along its edge, `S(s, r) = C(s) + r d`, with `C` one rational
quadratic. The corner is split about a centre `C` into one four-sided face at each place `Q_i` where
the contacts of blends `i` and `i + 1` cross, bounded by a half of either blend's end and by two
inner curves from the ends' splits `M_i` to `C`: the central split of an n-sided region [Charrot
and Gregory 1984; Chiyokura and Kimura 1983]. Each face is a rational Bézier patch of degree 6,
made directly; nothing is fitted.

#### Ends (derived)

Blend `i`'s end runs on its sheet at `r(s) = R(s)/W(s)`, `W` the section's
denominator and `R` the cubic with the crossings' rail places and slopes such that the end leaves
`Q_(i-1)` along `d_(i-1)` and arrives at `Q_i` along `d_(i+1)`. It is a rational cubic lying on the
sheet exactly, with an exact parameter curve there. Both directions and both sheets' tangent planes
at `Q_i` are those of the face the two blends share.

#### Across an end (derived)

If a patch's second row of homogeneous poles is its first moved along
`d` at the same weights, `H_j1 = H_j0 + w_j0 l_j (d, 0)`, then its derivative across that side is
`6 (sum_j B_j w_j0 l_j) / W d`, along the sheet's extrusion, so the patch is G1 to the sheet for any
offsets `l_j`. At `Q_i` the two rules give the same pole, because either end's tangent there is the
other side's direction. G2 to both sheets is impossible at `Q_i`: along the shared face the two
sheets bend differently there (the corner obstruction of [G2 at a corner](BRepPatchExt.md#g2-at-a-corner)).

#### Halves

Each half end is brought to the form whose two first and two last weights are 1: a
Möbius change of parameter (standard form), then a common cubic factor `(1, 2 - w_1, 2 - w_2, 1)` of
numerator and denominator, which gives a sextic. The end is split where the halves' speeds agree, so
the inner curves have unit weights and the halves join C1 at `M_i`.

#### Across an inner curve

Let `T` be the curve's tangent and `D_a`, `D_b` the derivatives across it of the patch it ends at `u = 1` and of
the one it ends at `v = 1`. The condition `D_a + D_b = 2 cos(2 pi / n) t T` is the G1 condition of a split of
`n` patches about one point [Loop 1994; Peters 2000]. It is linear in the poles. Its first two coefficients hold
by construction: the offsets next to `M_i` are `l_i (1 - cos(2 pi / n) / 6)`, with `l_i` at `M_i` (derived).

The remaining coefficients are solved for together with the free poles, which are kept fair by the least
thin-plate energy of the nets. That is one symmetric system for all three coordinates, with one more row that
holds `C` at a point `c`; it is solved by singular values, because the centre's equations are dependent. Its
solution is `x0 + g c`, with one `g` for every coordinate. Its last coefficient makes the inner curves' tangents
at `C` an affine image of the regular n-gon's.

#### Centre

`c` is where the surface bends least: the least bending energy `integral (k1^2 + k2^2) dA` over the patches
[Moreton and Séquin 1992]. The thin-plate energy of the nets stands in for it only where a patch is parametrised
evenly.

On uneven corners the patches differ widely in shape (a narrow face's crossing lies far from the corner), and
the nets' energy pulls the centre towards the long patches, below the corner's top, where the patch at a near
crossing has to bend back sharply (measured: at the optimum the corner bends 4 times as tightly as its blends on
a member where the nets' centre gives 130). Only `c` is chosen so; every other pole stays the nets' fairest for
it, so the conditions above hold for every `c`.

The energy is a sum of squares of the shape operator's entries in an orthonormal frame, `k1^2 + k2^2 = a^2 + 2
b^2 + c^2`, at the points of a composite Gauss rule (4 x 4 squares of 8 x 8 points per patch: the minimiser
moves by less than 0.7 % of the corner's size from that of 8 x 8 squares of 6 x 6 points on the uneven family).

It is minimised from the nets' fairest centre, which is closed form along `x0 + g c`, by quasi-Newton steps
(BFGS from the Gauss-Newton matrix [Dennis and Schnabel 1983, 9.2]), each halved until it lowers the energy,
until none lowers it by more than its rounding. The result therefore never bends more than the nets' centre
does, by this energy.

#### Orientation

At a crossing `Q_i` a patch is tangent to the face blends `i` and `i + 1` share, and
its sides leave along both blends' edges in that face; so its normal there is that face's, up to
sign, exactly. Each patch faces out of the body as that face does. (Its normal at its middle, against
the sum of the faces' normals, decided this before; on steep uneven corners that normal is nearly
perpendicular to the sum, and a sign error refused 24 valid corners as folded.)
`l_i` is a sixth of the chord from `M_i` to the splits' mean: the first leg of a straight curve to
the centre's likely place. This is a judgement; in measurements the peak curvature of the corner
depends strongly on `l_i`. Choosing it by the same fairness instead left a split with no positive
offset on uneven pentagons.

#### Measured

| Measure | Result |
| --- | --- |
| Seam gaps | up to 1e-14 of the size, on regular apexes and pockets of four to six faces (rho 0.3 and 0.7) |
| Seam angles | up to 1e-7 degrees |
| Bend of the corner against its blends (regular apexes) | at most 1.9 times as tight |
| Bending energy (regular apexes) | up to 7.6 % below the nets' centre's |
| Uneven apexes (apex moved up to half the base radius, faces from 12 to 114 degrees at the apex, heights 6 and 14, sizes 1e-3 to 1e3) | all 144 members build, without folds |
| Smallest ratio of `\|S_u x S_v\|` to its largest on a patch | 0.020 |

The corner bends more than 10 times as tightly as its blends on 22 of them, all on one pentagon whose faces span
12 to 38 degrees at the apex at height 14 (at most 34 times); on the others at most 7.6 times. That pentagon's
crossing on its wide face lies near the corner's top; neither the centre alone nor centre and first legs
together (measured by direct search) bring it below 10.

A template of tangents at `C` from the splits' actual directions, with the matching non-symmetric condition
`t_(i+1) + alpha_i t_(i-1) = beta_i t_i` (`alpha_i` constant along a curve, so the ends are split at speed ratio
`alpha_i`), was measured earlier and delivered fewer corners.

Choosing the centre by rules from the corner's points (along the line from the boundary's mean to the vertex, or
the point nearest the lines along which the inner curves start), or the nets' energy weighted by the patches'
side lengths or taken relative to their Coons nets, did not lower the peak curvature on every member.

Charrot and Gregory, *A pentagonal surface patch for computer aided geometric design*, CAGD 1 (1984).
Chiyokura and Kimura, *Design of solids with free-form surfaces*, SIGGRAPH (1983).
Loop, *Smooth spline surfaces over irregular meshes*, SIGGRAPH (1994).
Peters, *Patching Catmull-Clark meshes*, SIGGRAPH (2000).
Moreton and Séquin, *Functional optimization for fair surface design*, SIGGRAPH (1992).
Dennis and Schnabel, *Numerical Methods for Unconstrained Optimization and Nonlinear Equations*,
Prentice-Hall (1983).

### Experimental planar boundary curves on exact sections

For constant-radius sections on two regular surfaces, the equality of their offset
centres gives three equations in four support parameters. A section point in a plane
gives a fourth equation, with the transverse section parameter as a fifth unknown.
At rank four the common zero set is locally a curve, by the implicit function theorem.
Continuation in these five parameters permits an end curve to return to the same
contact. Contact and finite-carrier restrictions are additional event equations;
simultaneous events retain all their labels. At a smooth support partition, a new
support chart can continue the curve from the shared event and its incoming direction.
These constructions are derived.

`tracePlaneEnd` is an experimental evaluator, not yet called by solid attachment.
It evaluates the original support surfaces, uses their analytic offset Jacobian,
and corrects the section/plane residual directly. Centred differences guide the
section derivative; a predictor-deviation fraction of one eighth controls sampling.
Those are numerical choices, not existence or approximation certificates. The
output is corrected points on a seeded component, not a certified spline edge.
Component discovery, global exposure and installation of these components in the
body's faces remain open. Variable-radius laws, cliff contacts and section-transition windows
are outside this prototype.

| Generated boundary checks | Configurations | Result |
| --- | ---: | --- |
| Sphere/cylinder, two radii and two plane positions, circular and curvature-continuous laws, three scales and two poses | 48 | Both endpoints on one contact; independent plane/contact equations within 1e-9 of unit size, and the torus equation for circular sections |
| Same construction over two finite rational circular profiles, circular and curvature-continuous laws, three scales and two poses | 12 | Shared support event and final contact; both carrier/contact labels retained at their common endpoint |

The exposure of a circular sweep is a separate requirement: points lying inside
any other swept ball must be removed. Adsul, Machchhar and Sohoni, *Local and Global
Analysis of Parametric Solid Sweeps* (2013), section 3, Definitions 21, 23 and 24,
distinguish the contact set, its trim set and the exposed envelope
([paper](https://arxiv.org/html/1305.7351v1)). The numerical boundary tracer supplies
no global occupancy proof or closed-solid construction.

The separate `planeBoundaryCurve` fitter produces a shared cubic spline edge and
pcurves on a nonperiodic spline sheet and a plane. Local transverse charts permit
turns in either sheet coordinate. Hermite spans are refined until the stored
curves agree with both stored surfaces within the existing length tolerance.
The sheet certificate substitutes the polynomial pcurve into the Cox-de Boor
recurrence on every intersected knot cell. Outward-rounded Bernstein coefficients
enclose the homogeneous numerator of the discrepancy and a positive denominator.
This composition and quotient bound are derived; the basis recurrence and convex
hull property follow de Boor, *Spline Basics* (2001), sections 2 and 10
([chapter](https://ftp.cs.wisc.edu/Approx/cagdhand.pdf)). The certificate concerns
the stored sheet, not its discrepancy from the implicit section family.

`BoundaryGraph` assembles directed edge occurrences into retained cycles without
inferring incidence from contact labels. Distinct parameter-domain ports preserve
both uses of a seam even when their spatial vertices coincide. The current end
constructors import their original cycle without changing its edges or their
order. Replacing those constructors with traced contact intervals and support
events remains unfinished; no stop construction has been superseded yet.

| Boundary approximation check | Configurations | Result |
| --- | ---: | --- |
| Rational and polynomial spline representations, with and without an internal knot | 4 | Dense independent evaluations enclosed by the whole-span discrepancy bound |
| Cubic defect vanishing at both endpoints and the midpoint | 1 | Whole-span bound detects the interior defect |
| Parabolic same-contact section, three scales and two poses | 6 | Shared edge and both pcurves within 1e-7 length tolerance |
| Disjoint cycles, reversed seam occurrences and open/ambiguous incidence | 1 | Occurrences preserved; unresolved incidence identified |

Contact storage now consists of ordered retained intervals, each owning two
parameter/vertex events and its shared edge. The contact builder visits each
interval; continuity and tolerance collection visit all retained edges and
vertices. Existing endpoint constructors still populate one interval. They do
not yet consume the traced event graph.

The internal `completeCarrier` construction retains the solved support chart
while continuing a regular constant-distance branch to finite support
restrictions. Each terminal ball solves the offset-equality equations together
with the crossed support coordinate. Reparameterisation retains these stations
before fitting; it does not grow a new branch from the original edge. Closed
circuits and unresolved terminal events are not completed by this helper. It is
not yet called by the operation's end constructors.

| Retained-carrier check | Configurations | Result |
| --- | ---: | --- |
| Finite perpendicular spline planes, three scales and two poses | 6 | Both terminal contacts reach their exact restrictions; two disjoint contact intervals retain shared vertices and support pcurves within 1e-7 |
