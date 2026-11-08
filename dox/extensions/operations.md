# Operations in TKAlgoExt

Each entry gives the construction, the mathematics it rests on, the sources read for it, what is
not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule
chosen for the code, not derived; the constants are named in the source next to their use.
The whole-operation CPU and quality comparison with stock OCCT is in [benchmarks.md](benchmarks.md).

## BRepBlendExt and BRepFilletExt: rolling-ball blends

*Problem.* Two supports S1, S2 meet along an edge; for a radius r (constant, or a law r(s) along
the edge) find the envelope of balls of radius r that touch both. The ball centre c lies on both
offset surfaces, c in O1(r) ∩ O2(r) with O_i = S_i + r n_i (the sign follows convexity); the contact
points are p_i = c - r n_i(p_i), and the blend is swept by a section curve between p1 and p2 in the
plane normal to the centre curve. G1 holds at p_i by construction (the section is tangent to the
ball). The offset surface is regular where 1 - r k_i is non-zero for the principal curvatures k_i
of the support on the side the ball is on, so a radius reaching the smallest concave radius of
curvature is refused rather than built [Jaiswal et al.]. A section curvature-continuous with its
support needs only the normal curvature of the section to match that of the support in one direction
across the join, by the Linkage Curve Theorem [Salvi and Várady], so the section's third pole
is set by the support's curvature; that rule is a section *policy*, which a caller supplies (see [Section policy](../../README.md#section-policy)).

*Construction.* The centre curve is marched from a start point by Newton correction of the two
offset equations with curvature-controlled steps [Jaiswal et al., after Patrikalakis and Maekawa];
exact sections are fitted by a spline with knots at the data's continuity breaks and interpolation
at Greville abscissae (Schoenberg-Whitney) [de Boor]; the sheet is attached by trimming the supports
to the contact curves and the stopping faces to the ends. The conic section is the rational quadratic
Bézier whose shoulder weight is sin(w/2) rho/(1 - rho) for a wedge angle w, which is the circular arc at
rho = 1/2 [standard conic form]. Supports are extended past their ends by continuing the last
polynomial span, which keeps shape and parametrisation [Hu, as restated by Zang et al.].

*Sources read.* Jaiswal, Ghosal, Gurumoorthy, "Circular blending", DETC 2001 (theory, algorithm,
results; Hermann 1995 and Lukacs 1998 are cited second-hand there); Salvi and Várady, "Multi-sided
surfaces with curvature continuity" (2014; statement of the Linkage Curve Theorem, appendix proof not
read); Zang, Liu, Lai, GMP 2008, sections 1-3 (Hu's extension; Hu 2002 not obtained); de Boor 1972.

*Connected centre-length fitting.* At a defective rail knot in a constant-size
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

Sources: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided
Design and Manufacturing* (2009), section 2.1, equations 2.6--2.11; Peters,
*Geometric Continuity* (2001), sections 2--3. Centre fitting and geometric
acceptance remain sampled numerical checks. The generated polynomial and
rational-sheet family checks common terminal jets through order two, unchanged
window exteriors, rigid transforms and scales 1e-3 through 1e3.

*Arc length beyond a trimmed rail.* A rail can outlive the edge that supplied its parameter
range: the adjoining supports still carry the ball beyond an endpoint. For a nonperiodic multi-span
curve, `GCPnts_AbscissaPoint::Length` sums only the intersections with the adaptor's continuity
intervals (`computeLength`, composite case); an interval wholly beyond the edge consequently has
zero reported length. Expand a copy of the adaptor's parameter range before asking for arc
length: its original spline and rigid transform are retained, including knots hidden by the edge
trim. Integrating one supposed end span without exposing those hidden knots is insufficient.
The curve itself is neither fitted nor extended by a new geometric construction.

For the inverse, use the unchanged solver inside the edge's range. Outside, start with the local
speed estimate `distance / |C'(t)|` and double the parameter interval until its integrated length
brackets the requested distance. This bracketing argument is derived here: on a regular branch,
`d(length)/dt = |C'(t)| > 0`, so the length is strictly monotone; if a finite parameter realizes the
requested length, a finite number of doublings reaches it. Non-finite arithmetic and a step with no representable parameter progress do not supply a bracket
and are refused. Outside the trim, integration and inversion use the existing scale-resolved
`limit::solved` tolerance; no seam or construction bound changes. Integration and root solves remain numerical, without an interval certificate; no global existence claim is made past a
singularity or where the continued curve has insufficient length.

The closed-form check uses (t,t^2), whose primitive is
`t sqrt(1+4t^2)/2 + asinh(2t)/4`, across both endpoints and an inserted knot, both directions,
rotated and at scales 1e-3, 1 and 1e3 (294 interval/inverse checks). This primitive is derived here.
A held-out rational quadratic circular arc, also with an inserted knot and trimmed endpoints,
uses radius times the difference of endpoint `atan2` angles as its oracle (147 more interval/inverse
checks at those scales). A C1 curve that is a line up to t=1/4 and (t,(t-1/4)^2) beyond it tests
both trims that hide that knot, in both directions (98 interval/inverse checks); the earlier
single-span extension integration fails this family. No oracle uses the integration code for its
expected length.

This is the arc length of the **rail**, not a replacement of the centre path by an arc-length
march. Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and
Manufacturing*, hyperbook edition (2009), [section 5.8.2.3, equations (5.101)-(5.105)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html),
give the normalized intersection tangent and parameter ODEs; transverse intersection is required,
and the normalized cross product fails at tangency. Their topology discussion also requires a start
on every connected branch before marching. Jaiswal, Ghosal and Gurumoorthy,
[*An Efficient Algorithm for Circular Blending of Edges*](https://mecheng.iisc.ac.in/~asitava/sharad.pdf)
(2001), introduction and theoretical background, explicitly separate geometric spine construction
from trimming and integration into a valid body: that paper addresses the former. Neither result
by itself proves ownership of every retained trimmed region or that a coupled section satisfies
its prescribed law; those remain separate acceptance requirements.

### Fixed-contact conic windows

A final closed offset-circuit route can now attach a non-circular conic skin for a
cylindrical cavity meeting a transverse plane, at constant circular radius
`0 < r < R`. It follows every previous whole-request route. Neither previously
successful requests nor the crossing and solid acceptance checks are bypassed.
The carrier parameter frame is immutable: restricting a piece at a support seam
must restrict that circuit, rather than map a new complete circuit onto the piece.

Let `q` be the cylindrical contact's signed phase from the acute radial direction,
`A=R+r`, `h=2 asin(r/A)`, and `x=q/h`. The complete ellipse double-point domain
has half-width `asin(sqrt(((r/A)^2-c^2)/(1-c^2)))`, where `c` is the absolute
cosine between the plane normal and cylinder axis. It is enclosed by `h` when
`r/A >= c`. These formulas and the enclosure are derived; they are not a
published construction theorem. Outside `|x|<1`, the circular generator remains
unchanged. Inside, replace its shoulder `M` by
`M + (1-x*x)^2 (M0-M)`, with `M0` the intersection of the cylinder-axis line
through its cylindrical contact with the plane. Both shoulders lie in both
contact tangent planes. Positions and tangent planes at both contacts remain
exact analytically. The window and its first derivative vanish at the two ends,
so the retained/replacement first jets match. The existing homogeneous fitter
receives these sections and both contact pcurves together; knots of multiplicity
four in its degree-five longitudinal basis represent the actual C1 window breaks.
Homogeneous subdivision into two transverse spans follows the shoulder change.

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

After end cuts, mitres and ridge replacements restrict each continued contact to
its retained intervals, acceptance checks the support pcurve attachment actually
stores: joined where pieces share a contact, and projected from the exact sheet
where it replaces the fitted endpoint map. Acceptance intersects this UV curve
with every trim boundary and classifies each interval between events, including the contact endpoints. Ends
of tangential overlap segments are events too. On periodic carriers, conservative
UV boxes determine every translated boundary that can intersect the contact;
each witness belongs if it is in the face interior or within the face's own
model-space tolerance of a trim evaluated on the carrier. Boundary distances use
point/curve extrema and both curve endpoints, including degenerate pole edges.
The lifted trim is essential: the independently stored 3D edge can differ from
it by the edge tolerance. This union is the definition of membership in the
closed face's tolerance band (derived); it is invariant under periodic parameter
translations and does not turn a length tolerance into an angular tolerance
using a speed measured away from the witness. Membership is
constant between exact boundary events outside the tolerance band (derived from
continuity). The event locations and classification remain numerical OCCT operations, rather than an
interval certificate. Neither a carrier extension nor a missed or incomplete
intersection computation establishes finite support membership. A contact
outside its finite support is refused before attachment.

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

After the existing whole-request constructions refuse, a separate carrier may be obtained
by intersecting the constant-distance offset supports with `GeomAPI_IntSS`. Canonical offsets
are represented exactly through `Geom_OffsetSurface::Surface`; other offsets retain their
surface evaluator. The source-edge adaptor and its parameter range are kept separately.
The branch must be a closed circuit through the existing middle station, with compatible end
tangents. Its orientation follows that station. A closed clamped SSI spline is evaluated by
wrapping its finite domain, never by continuing its end polynomial. Every ball is still corrected
on the original offset equations. These numerical tests identify a branch, not a proof that all
branches were found. Open chains, carried ends, variable distances and cliff contacts are not
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

*Radius-law extrema.* Positivity and the scale used by the tolerances now come
from outward-rounded enclosures of every cubic span, rather than uniform samples.
Writing a span over x in [0,1], with h its parameter width, endpoints y0,y1 and
stored second derivatives m0,m1, its Bernstein controls are derived here as
b0=y0, b1=y0+(y1-y0-h^2(2m0+m1)/6)/3,
b2=y1-(y1-y0+h^2(m0+2m1)/6)/3, b3=y1.
The derivative controls are 3(b1-b0), 3(b2-b1), 3(b3-b2). Their interval signs
certify monotone subspans, whose extrema are endpoints. Otherwise de Casteljau
subdivision isolates every possible derivative root, including multiple roots,
until the coefficient intervals overlap or the parameter interval has no
representable midpoint. The remaining coefficient hull encloses its stationary
values. Thus every span endpoint and stationary region contributes to both bounds;
no fixed sampling width can miss an excursion. Nonfinite bounds or an unresolved
positive sign report a radius-law failure before any ball construction.

The range enclosure follows the nonnegative partition of unity of the Bernstein
basis; the derivative and subdivision argument above is derived here. Read:
Fuenfzig, Michelucci and Foufou, [*Polytope-Based Computation of Polynomial Ranges*](https://cfuenfzig.lima-city.de/Publications/PolytopeBasedComputationOfPolynomialRanges_CAGD2012.pdf)
(2012; preprint 2011), sections 3.1-3.2: interval extensions and Bernstein convex-hull
bounds. The existing interval primitives enclose arithmetic on stored binary64
inputs under round-to-nearest without fast-math. These bounds concern the stored
spline, not uncertainty in the input data or the linear solve. Native checks cover
a narrow negative natural-spline span, 18 quadratic extrema across independent
span widths and signs, three multiple derivative roots and a constant law.

*Steep and flat wedges.* In a wedge of opening psi between two planes the ball's centre lies
r / sin(psi / 2) from the edge (the centre is r from both planes, on the bisector), and its contacts
r cot(psi / 2) from it. Past an edge's end the blend is carried on while the ball goes on with the rail;
the carry stopped where the centre was more than 4 r from the rail, which a ball in a wedge sharper than
2 asin(1/4) = 29 degrees is already at the end, so its contacts stopped short of the faces that stop the
edge. The bound is now that many times the ball's own distance at the end over sqrt 2 where the faces
meet more sharply than at a right angle (unchanged at a right angle and wider), so it no longer depends
on the opening. Measured over straight edges between two planes, a plane and a cylinder, and two
cylinders, convex and concave, openings 5 to 175 degrees in steps of 5, rest, turned and scaled by 1e-3
and 1e3 (5040 blends, Round, Conic 0.7 and Chamfer; the plane-plane family used to find the cause, the
cylinder families held out): 516 refused before, at openings below 30 degrees, now 12, all of them
chamfers of setback 1 at a concave opening of 170 or 175 degrees beside a cylinder of radius 10, whose
chord between the feet runs through the cylinder (refused by name, as before). The planes' removed or
added volume equals r^2 (cot(psi / 2) - (pi - psi) / 2) per unit length (Round) and d^2 sin(psi) / 2
(Chamfer) to 3e-12 of the body.

Where the opening turns along the edge (the rim of a hole of radius R leaning a from upright through a
plate) the centres lie on an ellipse of semi-axes R + r and (R + r) / cos a, whose radius of curvature at
the acute end is (R + r) cos a. The rolling-ball surface (the envelope of the balls, a canal surface)
is regular where 1 - kappa r cos(alpha) is non-zero, alpha the angle of a point of the section from the
path's principal normal (derived here: along the path the surface c + r (cos alpha N + sin alpha B)
moves by |c'| (1 - kappa r cos alpha) along the path's tangent, its other terms lie along the section;
it is the regularity condition of a pipe surface in Maekawa, Patrikalakis, Sakkalis and Yu, CAGD 15,
1998, not read here); on the rim the arc between the
contacts passes through that normal at the acute end, so the surface crosses itself once
r >= (R + r) cos a: from a = 70.5 degrees for r = R / 2, 78.5 for r = R / 4. The contacts themselves do
not fold there (their speed along the path is s'^2 + R (R + r) |u'|^2 over |c'|, positive).
Just short of the crossing the fit missed its tolerance (lean 75 for r = R / 4): there the section
turns by up to 165 degrees, and a single rational quadratic's middle weight cos(turn / 2) falls to
0.13 with its apex 7.6 r off the arc, which the fit cannot follow as it moves. A conic section that
turns by more than 120 degrees anywhere along the edge (the faces' normals at its points; along a
whole chain) is therefore made of two spans, split at its middle by de Casteljau on the homogeneous
poles, the same curve (the 120 degrees, where the weight is a half, is a judgement; the 5040 straight
wedges are unchanged to 1e-9 of their volume). Where the surface does cross itself, the fitted sheet's
sections run back against the way the balls roll (S_u . (n0 x n1) < 0, sampled 65 x 7 places): refused
by name. Measured (the same 16 rims): the reference modeller makes the exact rolling ball
(8 faces) at 70 degrees and at 75 for the two rims with r / (R + r) < cos 75, and splits the blend
(11 faces, off the ball by up to 0.34 mm) for the other two at 75 and all at 80 and 85: exactly the
criterion above. The blend now makes all four at 70 and the two at 75 as the balls' envelope (every point r
from the centres' ellipse to 1e-6), and refuses the crossing ones; trimming them is open. The reference modeller
refuses G2 on every one of these rims.

*A size that changes along the edge* (`Request::radii`, `blendRange`). As the reference modeller makes it,
measured on its own results (its exported bodies of variable rounds, conics and G2 blends, sectioned across
the edge): the size r(t) is the cubic spline through the given places with no curvature at its ends
(the natural spline: the reference modeller's peak 1, 3, 1 on a box edge is r = 1 + 6 s - 8 s^3 on its first half
to the last digit), periodic round an edge that closes; and every section is the section of one size,
r(t), in the ball's own plane across the edge, the ball of that size resting on both faces with its
feet there (the reference modeller's rounds are circles there to 1e-12). The blend is G1 to both faces: at a foot
its tangent plane holds the section's tangent, normal to the face's normal in that plane, and the
contact's tangent, which lies in the face. Past the edge's ends the spline goes on along its end
tangents (its second derivative is zero there, so it stays twice differentiable), and the ball is
carried there only while its size stays within half and twice its size at the end (a judgement: a
steep law reaches nothing within twice the blend's size); the carry is twice the size at that end. The
places are breaks of the fit (the spline is C2 there and the stations follow its value). Not taken: the
rolling-ball envelope of balls of radius r(t), whose characteristic circle lies in the plane
(x - c).c' = -r r' (the derivative of |x - c|^2 = r^2; the feet lie on it since n.c' = r'); it was built
first and differs from the reference modeller by 0.06 mm^3 on a box edge with a linear law, and with the cubic of
`Law_Interpol` (OpenCascade's variable fillet's interpolation, the application's variable round's) by up to
2.4 mm^3 on a peak. A chamfer's setbacks follow the same law (the reference modeller has no variable chamfer: it
keeps the base distance or refuses); at an angle off a reference face the other setback is
d sin a / sin(theta + a), theta the opening between the faces' tangent planes at the edge's point (the
triangle of the edge's point and the feet). Distances off the two faces of a round, conic or G2 blend
must be equal where the size changes (refused by name otherwise); such edges add no tangent edges and
take the plain construction (no island, re-blend or mixed-convexity routes). A range keeps the change
the blend of the edge alone makes between the planes square to the edge at its limits, as
`filletRange` does for the round.
Where the premise below holds, the limited change is built directly (`directPremise`, `directRange`).
The premise: the blend is this module's own, of the one edge, whose two faces each change into one face;
the blend makes one new face (the sheet) that meets each of them in one edge (its contact); both limits
lie strictly inside the edge; every other edge the blend changed or made lies strictly outside the slab
between the limits, and each contact enters and leaves the slab once, through the two limits (sampled, as
the range's other checks are; not an interval certificate). Within the slab the change is then the strip
of each face between the edge and its contact and the sheet between the contacts, so each limit plane
meets it in three curves: across each face from the edge to the contact, and across the sheet between
the contacts. The result is the body in which, in each face's wire, the edge's piece within the slab
gives way to the curve across the face, the contact between the limits (the blend's curve and both of
its curves on the surfaces, restricted in range only) and the curve back; the sheet between the limits
enters with the side the blend gave it; and each limit is closed by a planar cap bounded by its three
curves, the cap using the sheet's section in the opposite sense. A plane meets a plane, cylinder, cone,
sphere or torus in lines and conics, taken exact (OpenCascade's analytic intersection) and their curves
on the face and in the cap by projection, which is exact for them. On other faces each sample of the
section is one Newton step sequence across, in the face's parameters, onto the plane; where the solved
points lie on an iso curve of the face to rounding, that iso curve is the section; otherwise the curve
on the face and the curve in the plane are interpolated at the same parameters, refined while they come
closer, and the edge takes their measured distance as its tolerance. Each new vertex takes the largest
distance from it to the ends of its curves and of their curves on the surfaces, at least 1e-7, and the
construction is refused where a gap exceeds the tolerances its curves carry (a construction that does
not close). Every other face, edge and vertex is the body's own, and the faces go back into the body's
shell. Outside the premise (the stock builder, a limit at an end of the edge, a blend of several new
faces), the change is isolated by Boolean operations as follows.
Where the blend reports which
faces it kept unchanged, those pairs cancel and only the remaining skins are intersected, once, in one
non-destructive intersection structure shared by both oriented differences (the separation of
intersection and building, Open CASCADE SAS, *Boolean Operations*, sections *Intersection Part* and
*Safe processing mode*). The retained pieces of the source faces are substituted in place and the new
sheets enter reversed; the pieces are grouped by their shared topological edges, never by a
tolerance merge or sewing, so the source surfaces and the edges' curve representations are kept as
they were. Each group is then a closed shell of the removed (positive signed volume) or the added
(negative) material. **Derived here:** that grouping gives whole outer shells when both bodies are one
solid bounded by one shell and share an unchanged face. A cavity shell of O\N with a face on the old
boundary would put points outside O inside a region enclosed by O\N, so O would have a void; one
with all faces on the new boundary would hold all of the connected body N, the shared face included,
which then lies on the old boundary inside that region: the same contradiction. N\O follows with the
roles exchanged. Bodies outside that domain, and the stock builder, intersect the whole bodies instead.
The stopping-plane clips and the final cut and fuse are Boolean operations.
Measured, direct against Boolean isolation on the generated family (box edge and cylinder rim, convex
and in a cavity, conic, G2 and chamfer, scales 0.01, 1 and 100, limits 0.17 and 0.63): the same face
kinds, faces within 5.2e-8 both ways and volumes within 4e-12 of each other, except where isolation was
wrong: a large G2 box in a cavity it returned unblended (the direct change is the range's share of the
whole blend's change to 3e-11, as at the smaller scale, the sections being the same all along) and two G2 rims in a cavity whose trim curves ended 1.2e-6 past their vertices. Vertex
tolerances fall to the real gaps: 1e-7 at scale 1 where isolation left up to 1.1e-4.

Measured: against the reference modeller's results (box edge and the flat-to-curved edge of a half cylinder, laws
linear, falling, peak and valley, Round and Conic rho 0.3) the blend's lie within 1e-12 to 1.1e-6 mm both ways;
round a closed rim the alternating law within 8e-6, but a law of three uneven places (1, 2, 1.5) 0.04 off,
the reference modeller's closed interpolation being neither the periodic nor the natural spline (open). Native
(`BRepBlendExt_Test`): rounds on a box edge (linear, steep, peak, valley) lie on their circles to 1e-6
and remove (1 - pi / 4) times the integral of r^2 to 1e-9 of the body, at rest, turned and at scale 1e-3
and 1e3; a cylinder's rim (periodic) and, held out, the ring of a sphere on a plane on their circles to
1e-6; chamfers of a growing setback, also at 30 degrees off a reference face and within a range, remove
the integral of their sections to 1e-9. G2 in the application (the audit's variable
and range rows): seams G1 to 0.01 degrees and G2 to 1e-4 / r.

*Why spatial stations stay.* `filletStations` (`BRepFilletAPI_MakeFillet` with a `Law_Interpol` law
through sizes at points in space, cut at limit planes) stays beside the section law, which cannot express
what it takes: sizes by length along a chain of tangent edges or round a closed contour (the section law
holds one law per edge), and limit planes that are not square to a straight edge (its range ends at
sections across the edge). Where both apply, on a box edge, measured both ways from tessellation points
to the other body's mesh (deflection 2e-3): a constant size between two limit planes gives the same body
(4.8e-12 mm, the same volume) at 0.55 to 0.8 of the CPU time; one limit plane 7.8e-5 mm, the same volume;
sizes 1, 2, 3 at three stations 0.046 mm and 0.53 mm^3 apart (`Law_Interpol` is not the natural spline);
sizes at four stations within two limits 4.5e-3 mm one way; a constant size round an oval rim of four
Bezier edges 1.6e-3 mm, at 1.5 times the CPU time. G1 and G2 of these pairs were not measured.

*Connected limited contours.* Plane intersections of each contour edge partition
its parameter range into intervals of constant halfspace membership. The existing
length tolerance shifts the classification plane by that distance; topology is
judged at the same tolerance as before. Intersections use OCCT's [GeomAPI_IntCS](https://dev.opencascade.org/doc/refman/html/class_geom_a_p_i___int_c_s.html)
(Open CASCADE 8.0.1, 2026; `Parameters` and `Segment`), which reports the curve parameters of intersection points. This remains a numerical
intersection solve, not an interval certificate. Count retained runs explicitly;
merge the first and last only on a closed contour. This connectivity argument is
derived here. Two transitions alone do not imply one run on an open contour: its
ends can both lie inside. Native negative controls use the parabolic rim y=10(2t-1)^2,
at three stopping heights 1, 5 and 9; each leaves two separated end intervals.

*Over a side narrower than the blend, and a pinched body.* Where the contact on a face lies past the
face's far edge along the whole edge, the face is used up: a round's ball rests on the far edge instead
(a cliff), any other section is cut by the face beyond the far edge (a notch). A contact that lies on
the far edge itself (a side exactly as high as the ball is wide) is taken the same way, as the ball
there rests on the face and its edge at once. Round a closed edge (a cylinder's rim over its side) the
notch's edge on the face beyond is a loop of its own, which takes the far edge's place, and the blend's
seam section ends on it. Measured: the top rim of a cylinder of radius 10 whose side is 0.5 to 2 high,
rounded or chamfered by 1 to 3 at least as high, keeps three faces and removes pi times the integral of
R^2 - x(z)^2 (the arc by its angle, the bevel in closed form) to 1e-9 of the body (round) and to the
length tolerance over the bevel (its edge on the bottom is interpolated); against the reference modeller's meshes of
the same rims (side 1 to 8, radius 1 and 3) every round and chamfer within 6.5e-7 mm, and its G2 at a
side as high as the ball within 1e-3 mm^3 of its volume. The reference modeller refuses G2 where the side is lower
than the ball; the blend makes it (the G2 section cut by the bottom face, its seams G2 to 1e-14 / r).
A hole tangent to a face beside the edge pinches the body: the hole meets the face in a line along
which four faces meet. No blend passes it and the blend is refused by name, as the reference modeller does; with a gap of
0.25 to 2 every round, chamfer and G2 lies within 1.3e-6 mm (G2 1e-4) of the reference modeller's mesh.

*Where a contact folds on any face.* A contact folds where the blend's size passes the radius of the
face's curve it follows (a setback r on a cylinder's wall beside a rim leaning a degrees off the axis
folds once r passes the rim's geodesic radius of curvature on the unrolled wall, about R / tan a): its
run back lies between two cusps, and the contact crosses itself at X. The crossing is found in the face's
parameters whatever its surface (formerly planes only), and the blend is cut at the two sections through
X. A round's loop is patched by its ball at X. A chamfer is the outer envelope of its sections (on the
cap, the contact trimmed at X: the distance set), so its two sheets meet in their crease, S(ra, va) =
S(rb, vb) on the one sheet, traced from X as the two sections close to the sheet's swallowtail point V;
from V to the other contact the two faces meet along one section, smoothly. No patch. Where the loop's
tail leaves the face over the blended edge itself, the face-crossing scan goes on past the loop, which
is cut away. A setback's side of the edge is read from both faces' normals at the edge's point, not at
its feet (a foot a setback round a curved face has its normal turned by the setback's angle, which where
the faces meet nearly flat or nearly closed turns the feet's direction over), and a chamfer contact's
speed is read along the edge, the normal of its section plane.
Measured: the slanted rims of a 4 or 8 thick plate, leaning 70 to 85 degrees, hole 1 or 2, setback 0.3
to 1, are chamfered in all 16 cases: every contact point lies at distance r from the rim within 8.6e-8 mm
and no blend point lies beyond r; the removed volume matches the union of the section triangles
(Monte Carlo, 2e5 points) within 2.5 standard errors (0.6 to 0.8%). The planar-cap chamfers of the
trimmed-outline corpus keep their distance-field volumes.

*Beyond the reference.* Three constructions go past what the reference modeller makes, and stay: a chamfer whose
distances vary along the edge (the reference modeller keeps the base distance or refuses), G2 where a narrow side is
lower than the ball (h < r; the reference modeller refuses it), and a chamfer larger than a round it runs over
(the reference modeller does not make it). Each is held to the same definitions and checks as the shared cases.

*Unproven.* Every marching limit (step counts, stall windows, Newton tolerances), the carry's stray bound
(4 r, or 4 / sqrt 2 times the ball's distance at the end in a sharper wedge), the radius/curvature
refusal margins, the choice between arcs and spline sections at corners, and the fit's pole and
stagnation budgets. There is no published bound for a G2 conic or quintic section as a function of the support's
curvature [searched; none found]; curvature-continuous requests therefore refuse in this tree until a caller
supplies a policy. Where a result is delivered it is measured again from outside (seam gaps, tangent angle,
curvature difference) and refused if it misses; that measurement is sampled, not a global certificate.

*Stock.* `BRepFilletAPI_MakeFillet` (ChFi3d). See the table below.

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

Jorg Peters, *Geometric Continuity* (2001), section 2, defines continuity by
regular reparameterization, rather than equality of arbitrary parameter speeds:
<https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf>.
**Derived here by the chain rule:** orient a common seam as
`A(0,s)=B(0,s)`, and write the first transverse reparameterization derivatives
as `alpha,beta` and the second as `gamma,delta`. Then
`B_t=alpha A_u+beta A_v` and
`B_tt=alpha^2 A_uu+2 alpha beta A_uv+beta^2 A_vv+gamma A_u+delta A_v`.
The mixed derivative follows by differentiating the first equation along s.
The parameter map must remain nonsingular; in these coordinates alpha is
nonzero, with its sign fixed by the chosen orientation. Contact and seam data
meeting at a corner must be derivatives of the same map. A quintic spline with
interior multiplicity at most three is C2, but this alone imposes neither those
support constraints nor regularity or fidelity to a prescribed section family.

Salvi, Varady and Rockwood, *Transfinite surface interpolation over specific
curvenet configurations* (2013), section 3.2, equation (3.1), uses parabolic
ribbons `P(s)+d T(s)+d^2 C(s)/2` for curvature continuity; sections 3.2 and 4
distinguish compatible corner data from rational correction terms:
<https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf>.
This supplies a construction framework, not a convergence or error theorem for
the current section interpolation. In particular a fit's positional error does
not bound its second derivatives.

**Derived negative control:** the graphs `(u,v,0)` and
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

## BRepPatchExt: boundary-first fill

*Problem.* Given a closed loop of edges, each optionally asked to meet a neighbouring face to G0, G1 or G2,
and optional interior curves and points, build one face. Candidates are tried in order: a plane; an exact
rational continuation of a torus or ball; a four-sided surface with the sides as its edges (a Coons patch whose
inner poles are moved to meet G1 and G2); for any loop, a fitted surface over the loop's mean plane (or over a
chart of its own, below), trimmed by the loop, with a net refined where the rim turns fast. The fit minimises a
weighted energy of squared second and third parameter derivatives subject to the constraints. The normal
equations are positive semi-definite for a connected net; they are solved by a skyline Cholesky factorisation
with a multiplier outer iteration. Chart coordinates over the loop are cotangent-weight (Dirichlet-energy)
coordinates, injective for a convex boundary [Floater03, Thm 4.1 (statement read, proof not)]; for a non-convex
boundary injectivity is not guaranteed [Karni, Gotsman and Gortler; Weber and Zorin], so the result is looked at
inside its loop before it is returned. Corner compatibility at a vertex where several patches meet follows the
vertex-enclosure condition [HermannPetersStrotman09].

### Measures and regularity

*Complete curvature normalization (derived).* A seam compares the two second fundamental forms in one orthonormal
tangent frame. For a symmetric form `[[a,c],[c,b]]` the spectral norm is `|a/2+b/2| + hypot(a/2-b/2,c)`. The
absolute mismatch is the norm of the difference; the relative mismatch is that norm divided by the larger norm of
the two complete forms (capped at one). Only two exactly zero forms have zero normalization: the saddle
`z=k*x*y` at the origin has norm `|k|` although both diagonal entries vanish. This solves the 2x2 symmetric
eigenproblem with the normal-curvature characterization [PMC09 §3.3-3.4, (3.26)-(3.43)]. It is a normalization,
not an approximation-error bound.

*Parameter units (derived).* Surface-jet rank is tested on independently normalized tangent columns. Their cross
product is the sine of their angle; zero-length columns or a sine at floating-point resolution (eight machine
epsilons, a roundoff policy, not a geometric tolerance) are unresolved. Coordinates use the cross-product dual
basis when the Gram determinant is poorly scaled. The dual identities follow from `x.(y cross (x cross y)) =
|x cross y|^2` for unit columns; the second fundamental form then acts on physical tangent vectors [PMC09 §3.3].
This removes absolute UV-unit cutoffs without asserting regularity of near-dependent columns.

*Sampled normal rotation (derived diagnostic).* Opposite normals at two separated observations trigger midpoint
subdivision, not a fold verdict: a regular cylinder can turn by more than a right angle between samples under a
nonuniform chart (exact rational circle: [OCCT8 GeomConvert]). Subdivision continues until consecutive normals
agree in sign, a sampled differential loses rank, or UV resolution is exhausted; the last two are reported as
unresolved regularity. This is not a certificate of local regularity everywhere or of global injectivity; the
separate Bernstein regularity certificate on chart fills remains required. Automatic candidate selection prefers
a graph when a four-sided surface has coarse opposed normals; this is a preference, not a fold verdict, and the
rotating candidate is checked with adaptive resolution if the preferred construction cannot deliver.

*Mean-plane degeneracy (derived).* The centered polygon-area sum is unresolved only below `8*n*epsilon*perimeter^2`
(n the number of vertices), not below a fixed share of squared perimeter, which would impose an aspect-ratio limit
on exact planar rectangles. Each centered radius is bounded by the perimeter, so bounding the cross products and
accumulated additions gives this conservative first-order roundoff policy. It is no lower bound on geometric area
or on the existence of a fill; the later geometric and tolerance checks remain required.

*Planar collars and graph directions (open construction; derived).* A planar support with unit normal `n` is a
height over a plane with unit normal `d` exactly when `n.d != 0`. Orthogonal projection restricted to its tangent
plane has singular values `1` and `|n.d|`, so its inverse condition number is `1/|n.d|` and the height gradient
has norm `sqrt(1/(n.d)^2 - 1)`. The `0.8` graph-selection threshold is therefore a conditioning policy, not an
existence condition. Three mutually perpendicular support normals have best common margin `1/sqrt(3)`. A simple
projected boundary and a global single-valued interior must be established separately.
For a polynomial tensor graph, affine reproduction gives exact planar boundary jets when every height coefficient
whose closed basis support meets a planar rim is set to that plane's affine height at the tensor Greville
abscissae: local support makes the value and first two derivatives planar at the rim, degree five with simple
knots gives C2. Incompatible planes must not prescribe the same coefficient; refinement separates disjoint
collars but cannot resolve inconsistent jets at a shared corner. The identities follow from [PMC09 §1.4.1,
(1.58)-(1.61)] by polynomial reproduction in both tensor directions: with `xi_i = (t_(i+1)+...+t_(i+p))/p`,
`p (xi_i-xi_(i-1))/(t_(i+p)-t_i) = 1`, so by (1.61) `(sum xi_i N_i,p)' = sum N_i,p-1 = 1`; the clamped left
endpoint fixes the constant, `sum xi_i N_i,p = u`, and partition of unity reproduces any affine height.
The construction is not enabled: in a three-direction planar-rim experiment (4 761 to 68 121 coefficients met the
rim checks) the interior-curvature ratios were about 15, 11 and 10 against the allowance 4. These are failed
candidates, not a proof of infeasibility; exact planar jets do not establish interior fairness, topology or
agreement with the intended shape.

### Guides, witnesses and trim checks

*Guides with positional tolerance (derived).* A guide is a positional constraint, not an exact endpoint two-jet,
so an endpoint tangent or normal-curvature mismatch cannot veto construction: `(x,y,z)=(.01*t,50,h*t)`, `0<=t<=1`,
lies within `h` of the planar square for every `h<=tolerance` whatever its endpoint angle. Construction and the
final measurements decide feasibility at the requested tolerance. This asserts no existence for arbitrary guides.
Endpoint jet information is kept only as a hint explaining a fit that already failed its result checks; it never
rejects a delivered fit. A rejected adjacent fit marks angular and complete-form curvature misses as well as
positional ones, and the between-row side checks call the same measurement and easing alternatives as
acceptance (derived). Finite refinement does not prove existence or convergence for arbitrary boundary jets.

*Owned witnesses near curved trims (derived).* A fixed-resolution trim polygon can misclassify points near a
curved rim. An interior-distance witness uses the parity of boundary arcs certified to avoid the query point:
positive-weight Bezier hulls enclose each arc in a convex tube about its chord [PMC09 §1.4.2, convex hull
property]; if the tube excludes the point, straightening the arc to its chord cannot cross the point and
preserves crossing parity [HormannAgathos01 §1-2, (1)-(3)]. Only uncertain arcs are subdivided; uncertainty is
kept once hull deviation reaches the floating-point enclosure guard. A parameter rectangle is a witness only if
its centre is inside the delivered trim and no boundary enclosure intersects it. This is not interval arithmetic;
unresolved conversions, offset approximations and near-boundary points give no witness.
Adjacent owned pcurves can represent their common vertex differently. A parity chain adds connecting segments only
where the gap fits the enclosure guard; for larger gaps ownership comes from `BRepClass_FaceClassifier` on the
delivered face, not an invented closing curve. Before ownership is propagated over a nonzero rectangle, every
boundary-arc tube and the whole coordinate box of every endpoint gap, enlarged only by the guard, is excluded;
unsupported pcurves cannot supply a cell certificate, and points classified on the rim stay unresolved. The
propagation is derived: an open connected rectangle disjoint from the boundary lies in one component of its
complement, so one authoritative interior point suffices. It assumes the face's own ownership semantics and
certifies no malformed wire. It changes neither the cells selected for surface regularity checks nor proves
injectivity of a fitted face.

*Interior trim checks (derived).* Requested points are measured against the delivered face, boundary included.
Positive-weight guide splines are subdivided until their control hull has a tolerance-ball cover about a measured
face witness: with `d` the witness distance and `r` the hull radius, `d+r` bounds every point by the triangle
inequality; for a convex planar trim the tolerance neighbourhood is convex, so every hull pole suffices [PMC09
§1.4.2]. For polynomial guide and surface pieces a paired cover avoids tolerance-sized balls: surface cells are
restricted to an in-trim UV rectangle and its diagonal, with the certificate partitioned at the guide's knots and
at every crossing of a surface knot by the UV chord; exact knot insertion and restriction give polynomial
intervals [PMC09 §1.4.3, (1.74)-(1.76)]. For bidegree `(p,q)` the diagonal poles are
`R_k = sum_(i+j=k) C(p,i) C(q,j) P_ij / C(p+q,k)` [PMC09 §1.3.6, (1.56)]. After degree elevation the largest
distance between paired guide and trace poles bounds the whole difference curve by the convex-hull property. The
bounds are exact-arithmetic arguments; floating-point conversions and pole calculations are not interval
certificates. A local projection only supplies an in-trim upper-distance witness; failing that, the global
trimmed-face distance query is used. For non-spline curves requiring approximation, the converter's reported
maximum error is reserved in the budget (not an independently validated bound). Unresolved coverage is refused,
never accepted from finite guide samples; coverage can be costly on curved faces at tight tolerances.

*Reusing projections (derived).* A projection onto a construction's carrier is not an owned witness on the delivered
trim: it can seed a UV search, but the surface value and trim ownership must be checked. Finite samples do not
bound a guide between samples; positive clearance from the trim establishes ownership of an enclosed UV path, not
the spatial distance from the guide to it. The final cover may be skipped only given an already certified
whole-guide cover on the same delivered surface and trim, or a proof that every later operation preserves it.

*Local extraction.* Knot insertion preserves the polynomial image in exact arithmetic [PMC09 §1.4.3]. Both
constructors of [OCCT8 GeomConvert_BSplineSurfaceToBezierSurface] copy the input surface before restriction and
knot insertion, so a bounded conversion per visited cell is not guaranteed cheaper than one whole conversion.
Ownership tests and local projections, not conversion, dominated the cover cost of a guided polynomial panel;
neither projection-seed reuse nor lazy bounded extraction established cost neutrality (measurements, not
complexity bounds).

### Energy assembly

*Assembly (derived).* A quadrature cell touches local poles `(a, b)` in order `p = (firstU + a) countV + firstV + b`;
these indices increase strictly, so an element's lower triangle lies in the global lower triangle. Each given
pole's coordinate (or its dot product with the height direction) is read once per fit and each element's
local-to-global indices once. The energy does not couple spatial coordinates, so visiting one coordinate at a time
keeps each right-hand-side entry's sequence of subtractions, and each matrix entry receives the same additions in
the same order (zero entries stay omitted from the skyline). Only indexing and repeated projections are removed;
quadrature, reductions, constraints and refinement decisions are unchanged. Four fixed corners determine the
affine nullspace, so the zero-energy affine plane is recovered from perturbed poles.

### Open boundaries

An open boundary asked for Quad sided, whose runs and gaps make four sides with corners, is filled by the
four-sided surface untrimmed: the poles of a gap side are free like the inner ones and nothing is asked along it
(a side with more poles at degree 5 than a fitted net may have is first brought down to degree 5 within a tenth
of the tolerance, as an arc is). If that quad is not a face to keep, the boundary is filled as any open one is.
Where the two rims face each other across two gaps, the fill must turn from the tangent plane of one to that of the
other, often more tightly than either rim bends (the end of a tube onto the end of a plate). The inside is then
also judged against the *bridge* between the rims: along each crossing, the Hermite curve that leaves and meets
each rim square to it in its tangent plane, at the speed of the chord, with the face's bend that way where G2
is asked. It is built from the rims' places, tangent planes and bends only, so it is pose- and scale-covariant,
and a bulge, ripple or spike that leaves the rims unchanged leaves it unchanged. The fill may bend four times as
tightly as the rim or the bridge within a tenth of the crossing's length of the same place (both factors
judgements). With no gap, or non-opposite gaps, the reference is the rim alone.

*A section as reference.* After a guided height fit fails, a planar section guide between two disconnected rims
supplies a reference candidate. Its plane meets each rim at an anchor; quintic Hermite segments extend the guide
to the anchors. At a smooth rim the extension's first and second derivatives come from the support jet,
differentiated subject to the section-plane constraint; an unsupported rim asks for position alone. The reference
is the single tensor spline `S(u,v) = R(u) + C(v) - R(u*)`, adding the two curves' poles and subtracting the
anchor, without a pole-count refusal.
*Derivation and scope.* This translation member of curve-network blending contains the section algebraically: if
the support is `R(u) + Q(v)` with its transverse two-jet independent of `u`, matching the section's two-jet to `Q`
at the anchor gives the same surface two-jet along the whole supported rim. A rail whose projection derivative on
the section-plane normal stays away from zero and a simple regular planar section give a regular injective
translation surface; these are sufficient conditions for this class, not an existence theorem for arbitrary
networks. Opposite rims with independent parameterisations are measured against their original curves. The
candidate must pass the delivered-rim, support, guide, validity, interior and sewing checks; a failed reference
keeps the refusal. The curve-network Boolean sum and its conversion to one tensor spline are in [Siggel18
§3.2.1-3.2.2, (10), (14)]; support-derived second-order ribbons and their continuity conditions in
[SalviVarady14 §3.3.2, Appendix]. The section-plane acceleration correction and this restricted route are derived
here. The interior curvature allowance is the unproven judgement below; a compatible guide can still be refused.

### G2 at a corner

*Corner compatibility (derived).* Take a common oriented tangent plane at a corner p. For a regular C2 fill F,
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

*Ease at a corner.* Where two sides asked for G2 meet at a corner whose faces bend differently (an arm of a round
meeting a strip whose edge bends across the arm's axis), no surface twice differentiable at the corner has both
bends there. It has the one bend that the two edges and the tangent planes along them leave it: four numbers (each
edge's normal curvature and the turn of its tangent plane along it) for the three of one bend, solved by least
squares, with each face's bend across its edge weighted a thousandth as much, so that it decides only what the
four leave open (for an edge that runs straight through the corner, the bend across becomes the mean of the
faces'). The bend each side asks across itself eases from that one at the corner to its own face's as
`1 - (10 t^3 - 15 t^4 + 6 t^5)` of the difference over the ease length: the quintic whose first two derivatives
vanish at both ends, so the bend asked stays twice differentiable along the side and agrees with the other side's
to second order at the corner. Within the ease a miss is allowed a share of no less than the ease itself (near a
corner of straight edges both the bend asked and the surface's are near flat). The fill is held to the eased
bends; the per-edge curvature it reports is against the faces themselves, so the corner's miss is reported and the
fill is not refused.
*Ease length (derived).* If the faces' bends are met on both sides beyond the ease, the fill's bend has to turn
across the corner by the whole difference at every distance d from it (on one edge it is that edge's curvature, on
the other the other face's), so it changes by that over an arc of about d; a net follows it only where d is no
shorter than its span. The ease is one span of the coarsest net (the loop's size over 6), at most half its side
(a judgement). Where no fill that eases so meets the faces beyond it (three rounds round a box corner at
tolerance 1e-6; unequal rounds and their planes round a spherical cut), the refusal names the corner and says
G1 can be had. The G2 analogue of the vertex-enclosure condition is stated open in [HermannPetersStrotman09]; this
is the elementary part of it that a fill can be held to. G2 is then met on the middle of each side only; coarse
nets that meet an eased loop can be lopsided; and the ease of one coarse span is a judgement bounded by
measurement, not an optimum. A flat strip meeting a bending one is held beyond its ease to the flat G2 limit,
which no eased fill meets, and is refused naming the corner (the same corner statement is made by a commercial
fill's manual: the corner's tangent plane and curvature "may not agree to the existing surfaces'").

### Spatial fit and chart

*Over a chart of its own.* A fill over the mean plane is a height over it, with slope `tan(acos |n.z|)` where the
rim's faces have normal n; the plane route refuses a loop where some `|n.z| < 0.17` or where the loop seen along z
crosses itself. Some loops are a height over no plane: where the face normals along the rim sweep a closed half
circle or more (a wing nose, a window round more than half a cylinder), every direction d is square to one of them
by the intermediate value theorem, so every height over every plane is vertical somewhere on its rim. The nearest
plane is read exactly: with the signs s of n.d fixed near the best direction of a search over the half sphere (a
Fibonacci lattice of 4096 directions in the loop's frame, so it turns with the loop), `max_d min s n.d` is the
distance from the origin to the convex hull of the signed normals [Boyd04 §8.1.1, support function], found by
Gilbert's iteration [Gilbert66]; every answer is a direction's own least, hence a bound from below. A loop refused
by the plane route whose margin, read at 65 places of every side asked G1 or G2 with their ends, is below the same
0.17 is filled over a chart (`overChart` in `BRepPatchExt_Patch.cxx`, chart in `BRepPatchExt_Chart.cxx`); any
other loop is filled over the plane, and a loop the chart does not fill is refused as over the plane.

The loop is first spanned by a disc of triangles: an even lattice in the unit disc, the loop on its rim at even
steps of its length (160 places), whose places in space minimise the sum of squared umbrella Laplacians at every
node of the rim and inside it, the rim held (discrete thin-plate fairing [Kobbelt97; BotschKobbelt04]). Along a
side to be tangent to a face, a ring of nodes outside the disc stands in the face at each piece's middle, less the
way the fill leaves the side (away from the face, read in the face's own parameters) times the height of an even
triangle; with the rim's umbrella near nothing, the nodes inside next to it stand where the face's tangent plane
goes on, the discrete form of the first-derivative condition that two fixed rings give for the biharmonic
equation [BotschKobbelt04]. Elsewhere the outer node is free (a natural end). The umbrella is the Laplacian of the
lattice's own parameters, so the disc is laid out again, twice, as an even lattice in the chart it lays flat to,
where an even triangle is near even on the surface and the umbrella near its cotangent Laplacian [Pinkall93], and
faired again; from the first of these on it is pinned to 33 places of every guide where it was nearest them (each
pin counting 100 umbrellas, a judgement). It is laid flat by Boundary First Flattening with the scale held at none
along the rim [SawhneyCrane17 §5-6.1; SpringbornSchroederPinkall08 App. E]: angles kept everywhere and lengths along
the rim, so that the thin-plate energy of the parameters is near the surface's bending energy, which it is not over
a plane where the surface stands steep. The flattening has no injectivity guarantee for a non-convex boundary
[Floater03 Thm 4.1; SawhneyCrane17 §7]; a triangle laid flat folded over refuses the chart.

A net free in space is then fitted over the chart as over the plane, each side at its places in the chart, its
tangent plane asked as `N.S_u = N.S_v = 0` and its bend as `n.S_nn = II(S_n, S_n)` about the surface of the pass
before (both conditions on derivatives in space, none on how the chart lies), and trimmed by the sides' curves in its
parameters. Such a net can meet every tangent plane by not leaving the rim at all (`S_n` running to nothing), which
is no parametrisation. In a chart that keeps the rim's lengths, a surface to which the chart is conformal leaves the
rim as fast as it runs along it, `|S_n| = |S_t| = 1`; along each side asked G1 or G2 the speed along the way the fill
leaves is held within 0.5 of that, asked only where it leaves that band (an inequality, not a wish; the band is a
judgement, since the chart is the disc's, not the fill's). The net starts as the coarsest of the ladder and spans
are halved where a place is missed, as for the Adjacent fit, up to the poles a net free in space is let have
(`freePoles`, 640); once a net meets the limits, the next finer is fitted and the coarser kept where they lie
within `settled` of each other (the stop below). The surface is kept only where `S_u x S_v` is shown away from
zero on every cell used: per cell it is a polynomial of degrees (2p-1, 2q-1) whose Bernstein coefficients come from
those of `S_u` and `S_v` [Farin02 §5.8]; by the convex hull property a cell whose coefficients all have a positive
part along one direction has no zero, and a cell that does not is cut in four, at most four times (the surface form
of the Bernstein Jacobian-sign test [Gravesen14]).

*Cost forecast (derived).* `netSpacing` predicts a chart's cost by replacing higher derivatives by those of the
osculating circle: a circle of curvature k has arc-length derivative magnitude `k^p` at order p+1. This is exact for
the circle. Curvature at a point of a general curve or surface bounds neither its higher derivatives nor those of
its chart, so the forecast and its 640-pole cutoff are computation policy: no approximation-error theorem, no
necessary pole count, no proof that a refused fill does not exist. [SandeManniSpeleers20 §1.2, (2)-(3), §6-7] give
upper error estimates in the actual Sobolev seminorm and, for mapped geometry, derivatives of the map; those are
L2/Ritz estimates, not pointwise curvature bounds, and cannot be inverted into a lower bound on a feasible pole
count. Replacing the surrogate requires measured constructions and the same final geometry checks.

*Spatial curvature condition (open).* An unconstrained spatial fit enforces one transverse curvature equation,
with five outer linearizations. G0 and G1 tolerance do not bound the omitted tangent and mixed second
fundamental-form components. A complete correction must check the actual nonlinear form and the delivered
parameter curves, including places between fitting samples, and establish regularity and interior acceptance.
Support-derived second-order ribbons give local jets [SalviVarady14 §3.3.2], not a convergence theorem for one
fitted spline. Globalized constrained iterations and feasibility restoration are in [NocedalWright06 §18.4,
Alg. 18.3]; their hypotheses and convergence must be verified for the assembled problem, and merely adding three
frozen curvature rows establishes neither. Open: a constrained feasibility solve with complete curvature conditions
and measured convergence, then refinement where the verified geometric residual requires it. The cap of a cone
standing 80 degrees or more against its circle asks for 529 to 576 poles at tolerance 1e-4, within the 640; its
coarse nets are not shown regular, and it is refused. The chart's own map is regular (a folded triangle refuses it);
what is not shown regular is `S_u x S_v` of the net fitted over it, which the chart does not decide.

*Chart constants not derived (judgements).* 160 places round the loop, 33 pins per guide, two relays, the pin's
100, the band's 0.5, four cuts of a cell, the 4096 directions (they only seed the exact step), and the 256 bins over
which the spans asked are read across the chart; the Favard constant is taken at its upper end, so that the spans
are, as far as the bound goes, enough.

### Stop and ladder

A net that meets every limit is not yet the fill asked for: that is the fairest surface of all that meet them,
which the nets of the ladder approach as they become finer. Once a net meets the limits the ladder continues until
two nets in a row lie within `settled = 1e-3` of the loop's size of each other, read from twice as many places as
the finer net has spans each way inside the loop to the nearest place of the other surface; the coarser of the
two that meets the limits is kept, and where the ladder or the allowed nets end first, the finest that meets them.
This is a Cauchy test: where the nets converge at least twice as fast as the ladder refines (ratio `rho <= 1/2`),
the kept net is within `d / (1 - rho) <= 2 d` of their limit, d the distance between the two. A share of the size
makes the stop the same at every size and pose. A net free in space and a height over the plane minimise
different energies and approach different surfaces; they are not one sequence, and once nets free in space meet
the limits no height is tried after them. `settled` is a judgement (about the height by which a tilt of 0.1 degree
moves a surface over half the loop's size). The cost is at least one net beyond the one kept; nets that a coarse net
met only roughly come out finer. Where nets keep differing by more than `settled` the whole ladder is run.
Open: the Cauchy bound needs the distances to fall, and they need not: on some loops (an eight-edge G0 opening
with nets 2.38e-3, 7.0e-4 and 7.5e-4 of the size apart) the ratio is above 1/2, the nets stand on a plateau below
`settled`, and the stop is a judgement, not a bound.

*A ladder that does not settle.* The distances between nets may not shrink at all (nets 3.8e-3, 3.1e-3, 1.3e-3,
2.9e-3, 1.2e-3 of the size apart, up to the finest). A condition is met once it is within its reach: the multipliers
stop at the first shifts that put every condition there, and the correction then moves the net onto the curves by
the least movement of its poles. What a net delivers is one surface of the band of all that meet the limits, at which
that iteration and correction stopped, not the least bending of the band. The band has an inside (the limits are
inequalities) and as nets become finer nothing draws their surfaces to one point of it. What is asked is the
fairest surface that meets the limits, the least of the bending energy `integral (k1^2 + k2^2) dA` over the band.
So where no two nets in a row come within `settled`, of the nets that met every limit the one of least bending is
kept (`bendingInside`, read at the places the inside of a face is looked at; `k^2 dA` has no unit, so the choice is
the same at every size and pose to rounding). A ladder that settles keeps what it kept. The energy is the measure
of fairness and the places are those of the inside check; nothing is newly tuned. Pose dependence is measured, not
removed: no change in whether a case is refused was seen over rigid-motion variants.

### Reading up to the corners

A fill is measured at the 129 places of Measure Continuity along each side and, between the first and last of them
and the side's ends, at places that halve their way to each end down to `cornerInset` (1e-6 of the side; a fraction,
so the same at every size; a judgement of how near is at the corner). G0 and G1 are read there as everywhere; the
bend is read there too, but not within an ease the corner's faces are shown to ask for (see `agreeAtCorners`),
where it is read from the first place on. Fills that pass the 129 places can leave their rim near the corner.

Where two sides that are only to be met meet at a corner p, a surface holding both is tangent at p to the plane of
their two tangents: each edge lies on it, so each tangent lies in its tangent plane, and two non-parallel tangents
span it. The plane may stand far from the surface around it, but it is asked at p only: along each edge the
surface's tangent plane need hold only that edge's tangent, and turns as fast as it will. Where the edges run on
within the G1 angle the two tangents span no plane and nothing is asked. Within the tolerance a surface need not
be tangent to it, and asking it of every such corner bends fills that were within before; so a correction onto the
sides holds, besides its even places, those of the measure's places near such a corner that the corrected surface
misses by more than the tolerance, and finer nets of the ladder are asked the plane at that corner.

### Size

Lengths are shares of the loop's size or of the tolerance, with one exception by design: where a fill goes over
into planes and other faces at once (or into a cone alone) its rim and its edges' parameter curves are made to a
tenth of `Precision::Confusion()`, OpenCascade's absolute least distance, so that its edges keep that as their
tolerance and STEP gives the joint back. A fitted net moved onto its curves ends its ladder when the bounded rim
correction leaves its own places further off than the tolerance (the correction's own residual, no constant).
The orientation of two collinear pieces of one straight side is decided by the sign of the `orient2d` determinant
only where it exceeds the bound on what rounding its products and their difference can do [Shewchuk97,
ccwerrboundA], a bound relative to the numbers themselves and so the same at every size. At a corner where both edges
run straight on rounds eased to one bend, everything the corner's four numbers ask is rounding; their share of
disagreement counts only above the rounding of the faces' own curvatures, `2^-40` of the larger principal curvature
(`jetRounding`, a judgement: a few thousand units in the last place, far below any curvature asked). The margin is
thin: rounding alone came to at most 9e-17 of the larger curvature against the floor 9.1e-13, but the least genuine
disagreement met, 6.7e-12, is only 7.4 times above it, so a genuine disagreement eight times smaller would be taken
for rounding and its turns dropped. A floor derived from the faces' own representation (each derivative off by up to
`(gamma_n + eps) sum |N_i^(a)| |N_j^(b)| |P_ij|` over the poles where the face stands, carried into the normal
curvature) was measured six to ten times above `2^-40` on flat strips and is not used: strips interpolated from posed
points are rounded by the solve's condition, not by the poles' storage, and no bound from their representation tells
that. A floor at the flat share of the G2 limit would be the measure's own and is the structural path once the
fold of the net over the plane without those turns is fixed. Open: a saved opening of near-flat spline faces
(curvature 6e-7 on a loop of 38) is filled only when the turns asked by its small but genuine disagreement are
kept; without them its fill folds over, an open defect of the net over the plane, not of the corner. At 1e-3 size
the corner of three rounds (tolerance below `Precision::Confusion()`) is made and one unequal-rounds corner is
refused at 1e3. The switch from the height graph to the spatial net at a support normal of 0.8 against the mean
plane changes a cone cap's area by 2.1 % between 36.8 and 36.95 degrees, and an analytic cone's cap changes between
two shapes from 70 to 74 degrees; a cap over a circle on a cone is refused as too steep from about 80 degrees.

### Ordering, chart covariance and symmetric requests

Permutation invariance of an annotated boundary differs from isometry covariance of its fill, `F(g X) = g F(X)`,
with supports, guides and continuity orders transformed together. If `g X = X` a covariant construction preserves
that symmetry; a unique minimizer supplies this only when its feasible space and complete objective are also
preserved (convergence of a penalized local solve does not prove uniqueness).
*Frames (derived).* A unique directed frame cannot be covariantly selected at every symmetric request: if `h X = X`
and `Q(h X) = h Q(X)`, an invertible frame gives `Q(X) = h Q(X)`, forcing `h` to be the identity. A chart defined up
to an orthogonal coordinate change works when all later stages respect that change. [DymLawrenceSiegel24 §2.2-3]
discuss further continuity obstructions and weighted frames.
The proxy chart is not such a complete construction: changing its seam rotates its circular boundary against a
regenerated triangular lattice, and a mesh whose connectivity respects cyclic shifts and reversal still cannot
resolve an ambiguous guide preimage (one of two equally close triangles can break a symmetry). A concentric-ring
mesh with four triangles per annular quadrilateral removes the arbitrary diagonal but is valid only while its
homothetic rings do not fold. An experimental guide initializer averages triangle-foot chart coordinates with
weights `area * exp(-(distanceSquared - minimumSquared)/meanArea)`; normalized weights are invariant under
isometries and uniform scale and the average transforms under affine chart changes (derived). The mean-area
bandwidth is an empirical choice, not an error bound; the average need not be a nearest preimage nor lie inside a
non-convex chart. Neither is a production caller. [SawhneyCrane17 §4.2, 5] give the flattening, not a guarantee for
these mesh or correspondence choices.
A fixed tensor spline space is generally not closed under chart reflections: reflecting `u^5 v^5` about an axis at
pi/8 gives `(u^2-v^2)^5/32`, outside bidegree (5,5). Coordinate reversal requires reversal-symmetric knot vectors,
which arbitrary adaptive knots do not ensure. Complete total-degree disk polynomial spaces, including complete real
Zernike levels, are closed under orthogonal chart changes [GreengardSerkh18 §2.3, 4]. A separate spline
approximation needs its own budget: approximating a reflection-symmetric surface within Hausdorff error e implies
mirror distance at most 2e by the triangle inequality (derived), with no derivative or regularity bound.

### Unproven

The interior acceptance gate (a fault when the interior bends more than 4 times the rim's bend, with a floor of
2/size) has no derivation; a bound of the form `K sqrt(2D)` in the energy gap D was tried and rejected. The bridge only
replaces the reference that the 4 multiplies; the 4 and the tenth are judgements. Against the bridge, defects that
leave the rims unchanged are caught less often where an elbow turns than against the rims alone (gentler bulges and
spikes pass more often). The flat G2 limit, the rules by which a ladder that has met nothing yet gives up
(`stillFalling`, `same`, `unruly`, which compare misses of unconverged nets and so depend on rounding and
therefore on pose), the span and budget constants of the net, the `settled` share and the 0.17 slope limit are
unproven.

*Stock.* `BRepFill_Filling` and `GeomPlate_BuildPlateSurface`: a plate fit that accepts G0, G1 and G2 point and
curve constraints and minimises a different energy.

## BRepCurveExt and BRepSurfaceExt

*Curves.* Offsets of planar curves with three corner treatments (arc, chord, extension to the crossing) are made on a
copy moved to the origin and scaled, so OCCT's absolute 1e-7 lengths do not change the result with the curve's size.
The outline of a body along a direction uses the edges and the turning curve `n.v = 0` of curved faces, with the
places where visibility can change found exactly. Tangent circles are found by Newton in a frame of the curves' own
extent. Also: projection, rebuild with a prescribed degree and count, refinement, and regions of coplanar curves.
A spiral is the helix (angle `2 pi t`, radius `R + p t`, height `h t / turns`).

*Surfaces.* Alignment of the control row beside a shared boundary on the line through the boundary pole and the
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

subject to the continuity conditions and fixed remaining poles. Trace normalisation (taken before constraint rows or
the prior are added) and the [0, 1] domain are design choices derived here, not a published surface-specific
weighting law. They make w dimensionless and keep rigid-motion and uniform-scale covariance for a fixed basis; the
pole metric is not invariant under change of basis or arbitrary UV reparameterisation.
For displacement unknowns d from the current iterate Pc, add `w a` to each free coordinate's diagonal and
`w a (P0 - Pc)` to its right-hand side; P0 is held through all nonlinear continuity iterations. Cost is O(N) time and
storage for N total poles, with no new poles and unchanged band; the term is applied after the energy cache so no
prior is cached. The condition/multiplier solve and final seam check still apply: the prior cannot authorise a seam
miss, and a surface that already meets the continuity is returned as it is.
*Derived bounds.* H is positive semidefinite (it integrates squared spline derivatives), so for `a > 0, w > 0`,
`H + w a I` is positive definite on the free coordinates and a nonempty affine constraint set has a unique
minimiser. For minimisers P1, P2 at `w1 < w2` on the same feasible set, adding their optimality inequalities gives
`(w2 - w1)(D(P2) - D(P1)) <= 0`, D the squared pole displacement: total displacement is nonincreasing, while
individual poles or curvature need not vary monotonically. Nonlinear G2/skew iterations change the linearised
feasible set, so this proves each fixed-row problem, not global convergence or monotonicity there. The final seam
certification is sampled, not a whole-domain Hausdorff proof.

### Collinear alignment: units and rational boundary rows

Collinear alignment uses exact degree elevation and the union of the two normalized boundary knot vectors with the
larger multiplicity at each knot, and edits rows directly (no degree-five fitting restriction). Knot insertion
preserves the parameterized surface [PMC09 §1.4.3, (1.74)-(1.76)]; `Geom_BSplineSurface::IncreaseDegree` preserves
the geometry and `InsertUKnot`/`InsertVKnot` with `Add=false` raise a knot to the requested multiplicity
[OCCT8 Geom_BSplineSurface].
*Derived:* for unit row directions a, b at angle theta, `|a+b| = 2 cos(theta/2)`. The symmetric fold check compares
this dimensionless quantity with 1e-9 (a conditioning choice retained from the unit-scale algorithm); physical
short-row degeneracy is checked separately. Scaling model lengths cannot change this angular verdict.
Rational supports need clamped boundary knots; rows carry weights as well as Euclidean poles. After exact refinement
the moving boundary and next row receive the reference's corresponding weights; deeper rows keep theirs. From the
rational surface definition [PMC09 §1.5, (1.91)] (derived): write the boundary as `S = A0/W0`; its cross derivative is
`c (A1 - S W1)/W0`, c the positive endpoint degree/knot factor. If the reference next-row weights are a constant
multiple of its boundary weights then `sum N_i w1_i P0_i = S W1`, so placing the moving next poles at
`P0_i + r (P0_i - P1_i)`, common `r > 0`, gives the opposite proportional derivative. This proves tangency for the
constant-ratio option, rational curved boundaries included. The proportional-weight condition is checked to 64
machine epsilons relative (an arithmetic allowance); a nonproportional reference weight row with that option is
unsupported. Without the option, pole lines align but neither polynomial nor rational variable row ratios
guarantee exact G1; skew changes the target lines and is not covered by the zero-skew argument.
At very large coordinates independently evaluated exact isocurves and surface points can differ by more than OCCT's
default absolute confusion tolerance. The face construction may then carry the requested positional tolerance onto
its new face, edges and vertices if the default metadata does not validate, and rechecks BRep validity and the
interior check. No poles, knots, weights or curves are fitted or moved by this step and no tolerance larger than the
caller's allowance is introduced; the option is off for other consumers of the face builder.

### Rectangular sheet extension

Extension requires a valid single-wire face whose boundary is its parameter rectangle. `GeomAlgoExt::parameterRectangle`
checks all boundary curves against the four sides: line endpoints, or the complete positive-weight spline pole hull,
which implies the property by the convex hull property applied to each constant coordinate [PMC09 §1.4.2] (derived).
Area agreement cannot establish the condition, since a narrow notch can have arbitrarily small area and finite
depth. Unrecognized boundary representations are not assumed rectangular.

### Coincidence of projected curves

Five common points and equal lengths do not imply equal curves. Plane projection removes an edge only after
establishing a common normalized spline basis (degree elevation, knot insertion, optionally reversal) with identical
positive weights. If every corresponding pole differs by at most the modelling tolerance, the rational
convex-combination formula bounds the displacement at every parameter by that tolerance, in both directions since the
common interval covers both edges, with outward-rounded pole-distance arithmetic [PMC09 §1.4.2] (derived). Edges whose
correspondence cannot be established are both retained. Degree elevation and knot insertion are OCCT's floating-point
implementations; the enclosure bounds the resulting common-basis representations, with the usual rounding limits.

*Unproven.* Work budgets and halving depths of the separation proofs, the crease threshold, interpolation grid sizes,
and the choice of pole counts ("twenty-one for any target"). The sampled checks are not Hausdorff certificates.

*Stock.* Offsets: `BRepOffsetAPI_MakeOffset`. Outline: `HLRBRep` (its documented limits on tangent and edge-on faces
are what the exact event search avoids). The rest has no direct stock counterpart.

## BRepOffsetExt, BRepSkinExt, BRepBooleanExt

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

*Tangent edges.* If two faces have the same normal n along an edge c and move by the same d, the point `c + d n` lies
on both offset surfaces with normal n on both, so they touch along the offset of the edge and its curves on the faces
stay. This is tested, not assumed: the two faces' points of each edge parameter after the move agree within the
tolerance, or within the edge's own tolerance where that is larger. Two faces tangent only to an angle theta come
apart by `|d| 2 sin(theta/2) <= |d| theta` more when both move by d. An input G1 only to theta (its edge tolerance
that wide) has its edge kept, and the result's edge carries the input's tolerance plus `|d| theta_max`. `theta_max` is
found, not bounded (judgement): 8 samples in each piece of the edge where its curves on the faces are polynomial and
cross no knot line of either surface, each local maximum refined by golden-section search. A bound would need the
faces' curvatures bounded along the edge, and the hull bound on curvature was three orders loose on blend sheets;
sampling found the larger maximum than a 4 097-sample reference on every nearly tangent edge tested, while 33 even
samples fell short by up to 36 per cent. The intended budget is that edge's carried tolerance and the input vertex's
plus the same spread (the endpoint gap limit below is not covered); the result is never claimed more accurate than
its input.

*Sharp edges.* The new edge is the crossing of the two offset surfaces, solved in the plane through c(t) normal to
the old edge: four equations `S1(u1, v1) = S2(u2, v2)`, `(S1 - c(t)).c'(t) = 0`, by Newton from the old edge's curves
on the faces, so it keeps the old parameter t. The curve and the curves on the faces are cubic interpolants of those
sections, refined locally (a midpoint is added where an interpolant is further than a tenth of the tolerance from the
section there); refinement stops at 16 halvings of the first spacing or 4096 samples (judgement: a curve a cubic
cannot follow there has a kink). An offset surface is one order less smooth than its basis (its normal takes a
derivative), so where the edge crosses a line along which a moved surface is less than C2, or the old edge (whose
tangent sets the section planes) is less than C3, the curve is interpolated apart on each side and joined C0 there;
the crossing is found among samples and by bisection. Newton steps that do not reduce the residual are halved; where
two surfaces only touch and the residual stalls, a point within half the tolerance of both is accepted (judgement:
the other half is left to the interpolation, a tenth, and the input's own edges).

*Vertices.* On a tangent edge's (or seam's) moved curve, crossed with the face around the vertex that meets it most
steeply among those that meet another at a sharp edge there (all tangent: the vertex moves along the common normal);
between two edges of the same two faces, in the section of one of them; else where three faces cross (the three whose
normals are furthest from a common plane). A vertex then not on every face around it is declined by this local
construction; that does not prove the topology must change, since approximate contact continuations can also leave the
other support.

*Certificate,* apart from the construction: a valid solid; points inside every face and along every edge at the moved
distance from their sources' unbounded surfaces, on the side moved to (foot point by Newton, Gauss-Newton where the
Hessian is not positive definite); faces clear of each other (below); volume grown by an outward move and shrunk by an
inward one; and, where every face moves, nothing nearer the body than the distance. A failure of any falls back to the
engine's routes. At a pole, where `S_u x S_v` vanishes against the derivatives' size, the certificate checks the
distance alone (judgement: below a millionth of `|S_u|^2 + |S_v|^2`), as the normal's sign is not defined there.

*Faces clear of each other.* The general self-interference check (`BOPAlgo_ArgumentAnalyzer`) on the exact faces, on
whatever the separation pre-pass (`BRepSurfaceExt_FaceSeparation.hxx`, pole hulls and Bernstein bounds) cannot prove
clear. The pre-pass does not read offset surfaces, so it runs on proxies: a copy whose offset faces are splines within
`eps` = a hundredth of the tolerance of them, each face's tolerance grown by `eps`. Two proxies proved further apart
than their tolerances have their exact faces apart too; that a proxy is one to one says nothing of its exact face, so
each offset face also goes to the check by itself unless the pre-pass leaves it in a part that is checked. The result
keeps its exact faces. Skipping the check by the reach [Federer59 Thm 4.8] (`|d| kappa_max < 1` per face, faces sharing
no vertex further apart than both moved faces reach) was dropped: a rigorous kappa bound per rational Bezier patch was
sound but three orders above the sampled curvature on blend sheets, and the reach leaves faces that share a vertex
unproven.

*Offsets that load again.* The binary BREP reader builds an offset surface with OpenCascade's continuity check
(`BinTools_SurfaceSet`), which refuses a basis that is C0 at a knot unless it finds it G1 there, so a body holding the
offset of a blend sheet whose parts are joined C0 in their parameters cannot be loaded from a saved file. Such a
face's offset is built unchecked for the construction and, once the face's extent is known, replaced by a spline over
its box within a tenth of the tolerance (`GeomConvert_ApproxSurface`), its edges' curves, vertices' points and
continuity records moved onto it and their tolerances grown by the error. An edge whose old curve is a line or a circle
keeps one: the offset of a line on a plane or cylinder is a line, of a circle round a surface of revolution's axis a
coaxial circle, with the parameter kept.
*Read back as written.* A body read from a saved file has the directions it was written with, bit for bit
(`gp_Dir::Restored` in the BinTools and GeomTools readers; [corrections.md](corrections.md#gp_dir-restored)). Open: the
X and Y directions of an axis placement (`gp_Ax2`, `gp_Ax3`, `gp_Ax22d`) are rebuilt from the main and X directions by
cross products when read and come back exact only where those products are; a reader that keeps orthonormal axes as
written is not made.

### Local focal witnesses for a sheet's normal sweep

`BRepOffsetExt_Focal` encloses non-periodic rational B-spline and Bezier supports directly from stored poles, positive
weights and knots, on each nonempty knot rectangle. The homogeneous first and second derivative control nets use the
derivative formula [PMC09 §1.4.2, (1.65)-(1.68)]; their active-pole hulls enclose their values by nonnegativity and
partition of unity [PMC09 §1.4.1, (1.58)-(1.61)]. Each hull is intersected with interval Cox-de Boor evaluation on the
current rectangle. Homogeneous poles, location matrices, knot differences, quotient derivatives and fundamental forms
all use outward-rounded interval arithmetic, with no prior conversion or knot insertion in ordinary floating point
[PMC09 §4.9.1].
Let I be the first fundamental form and II the second form against `S_u x S_v`; principal curvatures are the
generalized eigenvalues of (II, I) [PMC09 §3.4, (3.41)-(3.43)]. **Derived here:** the normal sweep `S + s n`,
`0 <= s <= d` (signed), has no focal point on a regular parameter rectangle precisely when `I - d II` is positive
definite, since its eigenvalues in an I-orthonormal frame are `1 - d k_i`. A negative diagonal entry or determinant
proves an eigenvalue has passed zero somewhere along the segment; Sylvester's criterion proves the positive case
without subtracting nearly equal curvature roots. The signed curvature convention is opposite to the hyperbook's
offset chapter; the focal condition is the same [PMC09 §11.3.1, (11.18)-(11.25)]. Face orientation reverses d. A
crossing witness must lie strictly inside the trimmed face; it refuses thickening with kind `Curvature`, naming that
source face, before copying or building.
A regular verdict proves only local regularity within the support's knot rectangles: not separation of distinct
normal segments, nor continuity across knots. Periodic or other supports, zero derivative denominators, singular
parameter regions, exact equality and exhausted work remain `Unknown`, and a loose interval never causes a refusal.
The construction runs for unknowns, so complete general self-intersection exclusion is open. The breadth-first
subdivision budgets (4096 rectangles, depth 24) are work choices, not geometric thresholds.

### Elementary focal distances

The offset Jacobian loses rank at a principal focal distance [PMC09 §11.3.1, (11.24)-(11.25); §11.3.2, (11.28)]; with
this implementation's normal convention the factor is `1 - d k`. For a sphere or cylinder of radius R the inward focal
distance is exactly R, and these two supports are classified by comparing the signed physical distance with the stored
radius, without the relative allowance used for fitted blends (a radius-10000 sphere moved inward by 9999.5 remains a
radius-0.5 sphere). For a cone of semi-angle alpha, `r(v) = R + v sin(alpha)` and the principal curvatures are `0` and
`-cos(alpha)/|r(v)|`; endpoint radii and a possible zero crossing give the full focal range over the trimmed v
interval without interior sampling. For a ring torus (`R > r`) the principal curvatures are `-1/r` and
`-cos(v)/(R + r cos(v))`: the inward crossing is at r everywhere; for outward d the other crossing is
`cos(v) <= -R/(r+d)`, and endpoints with the contained multiples of pi give the cosine extrema over the whole v
interval (curvatures with the same normal convention: [Sage10 principal_directions, torus example]). These
reductions are derived here from the parametrizations and the first and second fundamental forms; arithmetic and
trigonometric evaluation keep floating-point limits. Horn and spindle tori keep strict sampled classification.
Classification of general fitted supports stays sampled with the approximation allowance: not a proven
curvature-error enclosure. A general enclosure separating uncertainty from a strict crossing is open.

### Regularity of network surfaces

A lattice sized from the number of knots does not cover short knot intervals, so rejection witnesses include every
actual knot interval. Acceptance uses the outward-rounded Cox-de Boor derivative enclosure on every knot rectangle
implemented by `focalSweep` at distance zero. The derived identity `det(I) = |S_u x S_v|^2` lets the cross-product
enclosure establish regularity directly without dependency loss in `E*G-F*F`; a zero-free interval normal proves
local immersion and unresolved cells receive no such verdict. Basis and convex-hull bounds: [PMC09 §1.4.2-1.4.3].
Periodic finite covers use the stored extended flat knot sequence and cyclic pole copies, matching OCCT's BSplSLib
evaluation without fitting. This proves local immersion only; a global injectivity certificate is open (the argument
analyzer can flag valid periodic boundary seams and does not establish self-intersection freedom within one face).
To reduce interval dependency for rational surfaces, translate homogeneous coordinates by a constant local anchor,
`(X,W) -> (X-aW,W)`: this changes the surface to `S-a` and leaves every first and second derivative unchanged
(derived). All subtraction and multiplication are outward rounded. At distance zero only first derivatives are needed.
Local-immersion queries have no work-count cutoff; the finite subdivision precision limit reports unresolved, never
regular.

### A circular tube's parallel body after its hole closes

For a complete ring torus with core circle C of radius R and tube radius `r < R`, **derived here**: `T = C + B_r`, so
its outward distance-d parallel body is `T + B_d = C + B_(r+d)`, by associativity of the Minkowski sum and
`B_r + B_d = B_(r+d)`. With `a = r + d`, at radial coordinate `rho >= 0` its meridian satisfies
`(rho-R)^2 + z^2 <= a^2`, since the closest point on C has the same azimuth. When `a > R` this is the circle's disk
clipped at the rotation axis. Revolving the circular arc through `(R+a, 0)` between `(0, +/-sqrt((a-R)(a+R)))`,
together with its axis segment, constructs the boundary without the self-crossing continuation: no surface fitting or
new tolerance is involved, and the exact hole-closure equality `a = R` stays on the engine's route. It applies to a
valid positive-volume solid with one toroidal face, which maps to the new skin; validity, positive volume and the
clearance check still apply. The volume oracle is the meridian's first moment,
`V = 4 pi [ z0^3/3 + R (pi a^2/4 + R z0/2 + a^2 asin(R/a)/2) ]`, `z0 = sqrt((a-R)(a+R))`.
This is a trimmed parallel *body*, whereas thickening keeps the input sheet's normal sweep; trimming arbitrary spline
offsets requires tracing all self-intersection components [PMC09 §11.3.2, (11.28)-(11.29)].

### Offsets of moved faces: domains, seams, pipes

*Extending the approximation domain.* The kept-face offset approximates the exact normal offset on the UV rectangle of
the constructed face, including any extension needed to intersect adjacent moved supports. It supplies that rectangle
through `GeomAdaptor_Surface`, whose bounds select an evaluation domain, not `Geom_RectangularTrimmedSurface`, which
rejects bounds outside a bounded basis surface. Approximation tolerance, degree and segment limits are unchanged. The
normal offset is [PMC09 §11.3.1, (11.17)] and a spline as polynomial spans [PMC09 §1.4.2, (1.62)]. **Derived here:**
evaluating the end span beyond its knot interval continues the same parametric support used by the extended edge
curves; it neither changes the offset definition nor justifies extrapolating a fitted result outside its certified
rectangle. Approximation-error, topology, clearance and self-interference checks still apply to that rectangle.

*Both branches of a periodic seam.* A seam has two curves on its face but one 3D edge; they may differ by a period in
either surface parameter. When the end faces move, a finite seam pcurve can need the same end-span continuation as a
tangent edge. Both branches are fitted independently at the one-tenth-tolerance interpolation bound, and the modifier
gets the branch selected by that edge use's orientation; binding both curves by face identity alone would overwrite
one with the other. Branches needing no continuation keep their route. The representation contract is OCCT's
`BRep_Builder::UpdateEdge(E,C1,C2,F,Tol)` for a closed face and `BRepTools_Modifier`'s separate `NewCurve2d` calls for
the two uses. The continuation argument is **derived** from the parameter curves and (11.17). It does not prove that
arbitrary continuation is regular or non-self-intersecting; the constructed body is checked. Closed form for a
cylinder with finite spline seam pcurves: volume `pi (r+d)^2 (h+2d)`.

*Exact pipes.* For an exact pipe `P(s,t) = c(s) + r e(s,t)`, e the unit radial vector in the normal plane of the regular
spine, the normal is radial while the surface is regular; (11.17) gives `c(s) + (r +/- d) e(s,t)`, a pipe about the same
spine (**derived** from the parameterization and normal [PMC09 §11.6.2, (11.111)-(11.115)]). The local condition is that
the shifted radius times the maximum spine curvature is less than one; global separation also matters [PMC09
§11.6.3.4, Thm 11.6.1]. A fitted tensor spline through circular sections is not thereby certified to be this pipe, and
the kept-face operation has no retained spine/radius descriptor for a general spline support; inferring one from a
nominal size and replacing the support would change the input surface and is not done.

### Clearance coverage and rims

Clearance carries separate coverage and resolution states. Empty input, failed extrema, or exclusion of every region
cannot become successful infinite clearance. Distances are measured between whole trimmed regions and the original
skin, boundaries included; an interior lattice alone gives no lower bound. Near openings a region is measured clear,
or wholly inside the opening exclusion, or subdivided. Work exhaustion is unresolved, never a measurement.
*Exclusion bound (derived).* For a closed set A the triangle inequality gives `d(x,A) <= d(c,A) + |x-c|`; a containing
box supplies c and a radius covering the whole trimmed region, so `d(c,A) + radius < wall` excludes the region. The
radius is rounded upwards and the modelling confusion added. Trimming and distance computations keep OCCT's numerical
tolerances: a coverage check using those solvers, not interval certification of their extrema. A failed whole-region
solve can instead use the lower bound `d(x,skin) >= d(c,skin) - radius`, subdividing until it resolves the requested
clearance; a failed point solve or exhausted subdivision stays unresolved.
Distance queries use independent copies of each trimmed surface. Spatial edge curves inherited from an import can
deviate from the surface within their stored tolerance, and measuring them would measure representation error as wall
thickness; query curves are rebuilt from each face's pcurves and support with `GeomLib::BuildCurve3d`, whose maximum
deviation is returned (unresolved or excessive deviation fails the query). For query curve error at most e on either
shape, `|d(A,B) - d(A',B')| <= 2e` (derived by the triangle inequality); subtract it from the numerical distance and add
it to an opening exclusion's upper bound, with e = one tenth of the modelling confusion (a numerical estimate, not a
Bernstein-certified enclosure). A stored edge tolerance [OCCT8 Boolean Operations, "Interferences"] is not a justified
allowance to reduce the requested wall thickness.
*Walls and opening rims.* A shell's clearance domain is the inward descendants of the retained input faces; descendants
reached only from removed input faces close the opening and are rims, not offset walls [OCCT8
BRepOffsetAPI_MakeThickSolid, `MakeThickSolidByJoin`: the first shell holds the input faces except the removed ones,
the second is its offset]. This is tracked by input-face identity through construction, not by surface type,
coordinates, a distance threshold or a hash; a descendant shared with any retained input face is still measured;
opening descendants stay in the history; an empty retained domain still fails coverage. An opening in a rectangular
pocket can give an annular rim slightly beyond the removed face, its outer trim at the requested offset and its inner
trim touching the original pocket walls: independent whole-face distances give zero there, no positive lower bound
exists, and a uniform lattice can miss that zero. For plane separation h and in-plane distance r from the rectangular
opening, `d_open^2 = r^2 + h^2` while the adjacent trimmed wall can be at r; a spherical exclusion of radius t retains
distances down to `sqrt(t^2 - h^2)`, not t (derived by Pythagoras, `h < t`). Subdivision cannot repair a difference of
domains. These are numerical distance checks [OCCT8 BRepExtrema_DistShapeShape], not interval proofs of extrema.

### Neighbour-bounded planar push-pull

[ZouFeng19 §5.2, Alg. 3] build auxiliary volumes bounded by the neighbours' extended underlying surfaces and the
moved face's start and end supports; a prism is unsuitable for a slanted neighbour. The convex planar specialization
is derived here: write the body as an intersection of outward support inequalities `n_i.x <= c_i`. Translating the
selected faces by v changes only their constants to `c_i + n_i.v`, and the result is their new halfspace
intersection, exactly the neighbour-bounded volume update; active constraints may appear or disappear without a
topology rule. Moving constants leaves the recession cone unchanged, so a bounded original stays bounded whenever the
target is nonempty. Triples of independent support planes contain every vertex, giving the temporary box; no artificial
box face may survive, and each surviving face is tracked through Boolean history to its defining support. Both
Boolean differences between the original body and its reconstructed support intersection must contain no solid (this
rejects nonconvex input), and the final topology must pass the validity check. The construction is attempted only after
the sweep rejects a slanted neighbour. Curved supports and nonconvex neighbour-bounded volumes are open.

### Skins and sweeps

Frames along a path are rotation-minimising by the double reflection method [WangJuttler08], with a fourth-order error
in the step; stock `BRepOffsetAPI_MakePipeShell` twists a profile only where the path starts square to it. A skin's
cross-section lines are natural C2 cubics, whose minimum of `integral |s''|^2` over H2 interpolants is proved in
[utilities.md](utilities.md). Closed sections are made fit by cyclic seam alignment: the cyclic cross-correlation of two
sampled rings is maximised over start shifts and orientation [Dogan25] (the rotation proof by Kabsch/SVD read; the
version here fixes rotation). Draft of walls that `BRepOffsetAPI_DraftAngle` does not take (it handles planes,
cylinders and cones) uses the ruled isocline surface over the wall's neutral curve, whose every tangent plane makes the
draft angle with the pull direction.

*Gordon network.* Gordon [Gordon71] as presented by [Siggel18 §3.2.1-3.2.3, (8)-(17)]: the profile skin plus the guide
skin minus the intersection-grid tensor surface; also [PieglTiller97 ch. 10, p. 485]. Equation (8) is the
compatibility condition: each skin and the grid must use matching interpolation parameters and bases, which degree
elevation and knot insertion align without changing geometry. Compatibility and rational representability limits of
`GeomFill_Gordon` apply. A guided skin through three or more profiles first constructs the compatible Gordon network,
then applies the G1/G2 end correction below, avoiding assembly and factoring of a dense plate system for a network
already representable by spline interpolation. The plate construction remains when the network or its correction
cannot be represented or fails its checks. An accepted plate can have a different interior and bending energy, so
validation compares the delivered surfaces; exact curve interpolation alone does not establish equivalence. Coplanar
networks with planar supports keep the plate construction (a plane has zero geometric bending energy, and rebuilding
its unit network rectangle would change an accepted carrier and trims); the model tolerance defines the coplanarity
search and the plate checks decide delivery. For curved networks the accepted interior is compared with a matched
reference surface in both sampled directions.

*G1/G2 end correction (derived).* Two surfaces sharing an edge and their tangent planes along it meet G2 when their
normal curvatures agree in one direction across the edge (linkage curve theorem [PegnaWolter92]; for the geometric
interpretation and its tangent-plane-continuity assumption see [PMC09 §8.1.3.5]). The skin already holds the profile,
so G1 needs only the way across, `S_v`, turned into the support's tangent plane, and G2 only the normal part of
`S_vv` set to the support's second fundamental form at `S_v`; both moves are along the support's normal. The
correction is `h1(v) m1(u) + h2(v) m2(u)`, h1 and h2 the quintic Hermite functions of the end (`h1' = 1`, `h2'' = 1`
there, the rest zero) that vanish to the second derivative at the nearest place of the next profile. The other
profiles are kept and the skin stays C2; `m1`, `m2` are the moves less their straight interpolation between the guides,
so the guides are kept too. Its sampled field approximation is not a continuous error certificate. Where a guide leaves
the end profile out of the support's plane, or with another bend, the surface follows the guide. Nets whose skin carries
weights are outside this correction. Exact interpolation does not imply less bending energy (a counterexample raised
the integral from 0.9364 to 1.2673), so this candidate needs a family-wide geometry and fairness review before it
replaces accepted results.

*Pairing of unlike closed sections (open).* Approximating one conic by non-rational cubic pieces is not the problem of
correspondence between unlike sections. [Floater97 §2, Thm 3.1] give geometric Hermite interpolants and a Hausdorff
bound for positive normalized middle weights up to 3 (for degree three the subdivision order is six); it controls the
conversion of each section, not which points on two sections are joined. Repeated bisection of three initial parameter
regions can explain a knot pattern without determining the correspondence. Rational reparameterization supplies none
either: [HuWangWang18 §2.1, (3)] change weights by `w_i -> lambda^i w_i` under `t = lambda u / (1 - u + lambda u)`,
`lambda > 0`, preserving degree and geometry; their §2.2, (7), (9), (12) use curvature maxima as breakpoints,
normalized arc lengths as new breakpoints and matching derivatives, as choices for one curve's speed, not a theorem
selecting points between unlike curves; a three-region knot pattern does not establish Möbius pieces (a one-map
diagnostic through three places per region missed the observed conics by 0.234 to 0.550 model units: evidence
against that hypothesis, not a bound on piecewise reparameterizations). Compatibility is distinct from choosing
correspondence: [PieglTiller97 §10.3, pp. 458-459, (10.7)] assume common section degree and knots and interpolate the
control points (homogeneous for rational sections); [EngleitnerJuttler18 §2-4] permit incompatible section spaces
through patchwork splines, with existence in the chosen space (Thm 3, under shadow-compatibility and
smoothness-monotonicity) and a constrained energy minimization (§4.3). In every case the sections'
parameterizations are given, so none identifies an unknown pairing from an output knot pattern. No pairing law is
installed.

*Solid publication.* Planar end caps are the only cap construction of the smooth section skin. A failed cap is a
construction failure for a solid request; an open sheet request returns its sheet. Before publication closure is
checked with `BRepCheck_Shell::Closed`, orientation with `BRepLib::OrientClosedSolid` and validity with
`BRepCheck_Analyzer` [OCCT8 BRepCheck_Shell, BRepLib]; disconnected shells are not silently discarded. These establish
topological closure and validity, not freedom from geometric self-intersection. Nonplanar end caps are open.

### Booleans

Cell selection from one General Fuse with a ladder of fallbacks, each rung accepted only if checked; a sheet tool
divides rather than subtracts.

*Unproven.* Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean
fallbacks, the sampling counts of the offset and skin validity checks.
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

Sources: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided
Design and Manufacturing* (2009), section 3.3, equations (3.26)-(3.28),
[normal curvature](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node29.html),
and section 8.1.3.5,
[linkage theorem](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node156.html);
Peters, *Geometric Continuity* (2001), sections 2-3,
[author manuscript](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf);
Salvi, Varady and Rockwood, *Transfinite surface interpolation over specific
curve-net configurations* (2013), sections 3.2-3.3, equation (3.1),
[parabolic ribbons and corner compatibility](https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf).
These sources establish the geometric continuity and ribbon compatibility
framework; the coefficient reduction and proposed correction above are derived.

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

A local diagnostic construction can share coefficient value, first derivative
and second derivative at every rail station. On a span of length h the quintic
Bezier coefficients are P0, P0+h D0/5, P0+2 h D0/5+h*h A0/20,
P1-2 h D1/5+h*h A1/20, P1-h D1/5, P1. Shared jets imply C2 in exact arithmetic,
even when redundant knots are stored. Floating-point cancellation on small spans
still needs a physical curvature check. For example, the endpoint second
derivative is 20*(B2-2*B1+B0)/(h*h). If each stored coefficient component has
absolute error at most e, its second-derivative error is at most 80*e/(h*h),
before evaluator rounding; two independent sides can differ by 160*e/(h*h).
These are derived componentwise bounds, not curvature bounds: the surface metric
and normal also enter curvature. A genuinely C2 knot representation avoids
independent endpoint storage, but cannot supply missing support constraints.
In a common polynomial tensor Bezier
basis, the maximum norm of corresponding control-point differences bounds
same-parameter displacement by the convex-hull property. It bounds neither
normal differences nor second fundamental forms. Knot insertion rounding must
also be enclosed before calling this a numerical certificate.

After changing a sheet, build its contact curves and stopping rims from that
sheet. A projection fit's error in support parameters is not a world-space
error bound when the support metric is poorly scaled. Check the composed
support curve against the actual sheet in physical units, including endpoints;
sampled checks remain sampled evidence. Exact one-sided jets where rail knots
reach a retained contact can reveal errors missed by interior edge samples.
Preserving a positional proposal allowance cannot replace the original contact,
section, topology, regularity and internal-curvature acceptance tests.

Sources: Hartmann, *Geometry and Algorithms for Computer Aided Design*
(October 2003), sections 11.3.2 and 11.4, pages 133-135,
[common offset-intersection parametrisation](https://www2.mathematik.tu-darmstadt.de/~ehartmann/cdgen.pdf);
Salvi, Varady and Rockwood (2013), section 3.2, equations (3.1) and (3.3),
[parabolic support ribbons](https://3dgeo.iit.bme.hu/papers/multisided/transfinite-g2.pdf);
Patrikalakis, Maekawa and Cho (2009), section 1.4.3,
[knot removal and continuity](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node18.html).
These sources supply constructions under their stated regularity and
compatibility conditions; they do not prove existence of an arbitrary
support-constrained approximation in a fixed tensor space. The trace identities,
Hermite diagnostic and preservation qualifications above are derived here.

#### Normal-acceleration corrections and numerical evidence

For a regular seam whose tangent planes agree, let t be its unit tangent and
c = n cross t. A correction delta*n to S_uu changes II(c,c) by
`delta/(S_u dot c)^2`, with the seam position and first derivatives held fixed.
This identity is derived from the inverse surface metric. Matching this scalar
alone is sufficient only if the along-seam and mixed second fundamental forms
already agree. First-order coupling of independently fitted pieces does not
establish those hypotheses; check the full form in one physical tangent frame.

A derived compact acceleration function on a one-sided interval of width h is
`h*h*x*x*x*(1-x)*(1-x)/2`, where x runs from zero at its outer boundary to one
at the seam. It vanishes through order two at the outer boundary, has zero value
and first derivative at the seam, and has second derivative one there.
Its maximum is `54*h*h/3125`. Thus a coefficient field bounded by Q displaces
the surface by at most `54*h*h*Q/3125` in exact arithmetic. This is a positional
proposal bound, not a curvature or support-contact certificate. Transverse
collocation of the field does not establish an identity between samples.
Successive knot insertion and world-coordinate coefficient rounding can destroy
the intended second derivative on very short spans. Adding more transverse
sites cannot repair that independent longitudinal conditioning problem.

Use independent, well-conditioned derivatives before judging a construction.
For a degree-p spline with expanded knots U, its derivative coefficients are
`p*(P[i+1]-P[i])/(U[i+p+1]-U[i+1])`, on the knot vector with the first and last
knots removed. Apply the formula in each tensor direction and to homogeneous
coordinates for rational surfaces, then apply the quotient rule. Translating
positions by a fixed origin before evaluation leaves all physical derivatives
unchanged in exact arithmetic. Form first homogeneous differences directly as
`w_b*(P_b-P_a)+(P_a-O)*(w_b-w_a)` for a fixed origin O (derived identity);
subtracting the origin from each pole before differencing can itself discard
the small difference on a short span. Higher differences then act on derivative
coefficients. Summing differentiated basis functions against large world positions
can instead cancel terms proportional to the inverse square of a short span.
Derivative control nets reduce evaluator cancellation; they do not recover
information already lost in stored coefficients. Increased-precision evaluation
of the same stored net distinguishes these two sources of error.

A positional fit across an entire carried rail is likewise insufficient: making
all internal knots structurally C2 can approximate a tangential parameter kink
using very small spans while introducing severe support curvature. A structural
alternative is to solve compatible contact ribbons and their parameter freedoms
jointly, controlling the required derivative norms and the longitudinal
conditioning before fixing a tensor partition. This is a design direction,
not an existence claim for arbitrary boundary data.

Sources: Nowacki, Kaklis and Weber, *Curve mesh fairing and GC2 surface
interpolation* (1992), sections 3.2-3.3, equations (26)-(42),
[quintic Hermite interpolation and compatibility](https://www.numdam.org/item/M2AN_1992__26_1_113_0.pdf);
Peters, *Geometric Continuity* (2001), sections 2-3,
[reparameterisation conditions](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf);
Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design
and Manufacturing* (2009), section 1.4.2, equations (1.65)-(1.68),
[derivative control coefficients](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html),
and section 1.4.4,
[tensor-product surfaces](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node19.html).
The compact correction, its maximum and the conditioning conclusions above
are derived here; the cited sources supply the interpolation and continuity
framework, not a proof of convergence of this correction.

### Face-offset endpoint certificate limits

A normal-offset surface [PMC09 §11.3.1, (11.17)] and a valid BREP do not by themselves give a bounded face-offset
result: the boundary must also realize the required incidences within an independently assigned budget. For every
oriented endpoint occurrence measure `|V - C(t_V)|`, and compare the 3D edge curves with both lifted pcurves. A
metadata tolerance enlarged afterwards to cover a measured gap is not an independent result bound. Closed edges need
both endpoint orientations. The certificate checks every oriented vertex/edge occurrence against
`max(operation tolerance, input vertex tolerance, carried input vertex tolerance + incident normal spread)`, assigned
before any incident curve is fitted; both lifted seam pcurves are checked after rebuilding and optional surface
approximation. The inherited normal-spread estimate is sampled, not an interval enclosure.

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

*Why an old trimmed-contact bound does not certify a continued contact.* Source lifts continued independently (by
homogeneous de Casteljau continuation of the end profile span, then rotation, which matches the offset `S + d n` to
floating-point error) already differ at the first refused contact by 3e-5 with normals differing by 1.6e-3 although
they agree at the old vertex. The triangle inequality gives `gap_after <= gap_before + |d| |n1-n2|` (derived). It
neither authorizes a larger output tolerance nor shows that a jointly constructed continuation cannot exist.
Closest-point stationarity: [PMC09 §7.2.1, (7.8)-(7.9)]; SVD least squares and Gauss-Newton: [NocedalWright06
§10.2-10.3, (10.20)-(10.26)]. Neither establishes that a locally refined vertex connects to the desired intersection
branch, a failed iteration is no impossibility proof, and reconstructing nearby exact section intersections and
minimizing contact distances in model space also failed the measured continuation; none shows that no offset exists
within the input uncertainty.
*Structural alternative (proposed, not implemented).* Construct each junction together with its incident contact
curves, track intersection branches before trimming, and carry explicit endpoint and edge/support budgets through the
fit. [Mukundan04 §2.2, (4)-(6), §3-4, Thm 1] distinguish isolated tangential contacts, branching and higher-order
contacts, then validate tracing and map parameter enclosures to model-space enclosures; their assumptions (an identified
branch, a bounded starting point, regular rational interval surfaces) make this no ready theorem for arbitrary
normal-offset surfaces or approximate contact seams. The [OCCT8 BRepOffsetAPI_MakeOffsetShape `PerformByJoin`] contract
distinguishes intersection joins from added arc faces and promises no self-intersection elimination. A general
construction must establish its own existence conditions and final bounds.

### Separated tangent junctions

The fixed-incidence construction runs first and keeps its successful results. After a vertex fails the
point-on-every-support check, a single-shell solid with three incident faces and a retained tangent edge can instead
split the junction. The tangent edge ends at its original endpoint transported to the offset; each of its two
continuations intersects the third offset support independently, and the two sharp edges use those separate endpoints.
The old support surfaces stay exact. A ruled strip `R(u,v) = (1-v) C0(u) + v C1(u)`, boundary curves on [0,1] [PMC09
§9.7.1, (9.74)-(9.76)], fills the boundary between the continuations and its transverse end closes the third face;
face lineage gives the strip all three incident parents. `BRepFill::Face` builds it through `GeomFill_Generator`. It is
a boundary closure, not a claim that the strip is a regular developable surface or an exact normal-distance envelope
(normal offset: [PMC09 §11.3.1, (11.17)]; separately generated edge and corner patches in solid offsets are discussed in
§11.1.4, which does not prove this strip's validity). Unsupported arrangements refuse.
*Checks.* The seam's two transported endpoints must agree within the operation tolerance; rail fits use one tenth of it;
the transverse chord must project to the third support within it. Every oriented endpoint and each new boundary's lifted
pcurve is checked against that bound independently of output metadata, and validity, self-interference, signed-distance,
volume and whole-body clearance checks still run. Curve/support checks are adaptive fits followed by sampled validation,
not interval certificates of the whole strip.
*Shell orientation.* Shell reconstruction keeps the orientation of each source-face occurrence, including faces outside
the split junction: the modifier's shape map identifies shapes independently of occurrence orientation, and its stored
face is not a replacement for the oriented use [OCCT8 Modeling Data, "Orientation"; `BRepTools_Modifier::Rebuild`,
"copy the sub-elements"]. The stored orientation of an unchanged face can give a shared edge two identical uses, so a
closed collection of individually valid faces is not an oriented boundary. Carriers, curve parameters and contact bounds
are unchanged.
**Derived diagnostic:** for two regular contact continuations with a common position and tangent, the normal separation
is `0.5 (k0-k1) s^2 + o(s^2)` in the common tangent arc coordinate s. The curvature is taken on the continued side of a
boundary knot (a symmetric finite difference across a second-derivative break averages the two sides). This Taylor
estimate diagnoses the failed old incidence and never assigns a result tolerance; it does not certify existence,
regularity or nonintersection of a finite strip.

### Transported contact uncertainty is confined to its measured domain

Let `p(t), q(t)` be corresponding lifts of a source edge and `n(t), m(t)` their outward unit normals. From the normal-offset
definition [PMC09 §11.3.1, (11.17)] the triangle inequality gives the **derived** bound for equal signed moves d:

`|(p+d n) - (q+d m)| <= |p-q| + |d| |n-m|`.

For unequal signed moves replace the second term by `|d1 n - d2 m|`; if the representative 3D curve rather than the other
lift is used, its discrepancy is added. A bound on the original closed parameter interval does not bound an
extrapolation outside it, so a large continued residual cannot justify itself by replacing the source uncertainty; and a
measured nonzero input error is no necessary lower limit on output accuracy (the inequality is an upper bound, not a
nonexistence theorem). The operation's 1e-6 allowance is not replaced by a constant near 1e-5: on 47 generated sources
with approximate tangent contact against revolution supports (input edge tolerances 1e-7) the sample maxima were 5.2e-8
for the source lift disagreement and 5.5e-8 for the pointwise transported bound at `d = 0.05`; sample maxima are
diagnostics, not certified suprema, and the endpoint and paired-pcurve checks remain necessary.
A reference result preserving the whole face-edge-vertex incidence graph, or returned without error, is no strict
geometric oracle: bodies were seen with curve endpoints within 7.6e-10 whose junction points were 2.3e-5 from the
exact offset revolution support, and others with endpoint discrepancies above 5e-5.
*Offset of a revolution by its meridian (derived).* Write the generating curve as `(r(v), h(v))` about the axis, `r > 0`,
`L = sqrt(r'^2 + h'^2) > 0`. For surface-normal offset d the offset meridian is `(r + d h'/L, h - d r'/L)`, by
evaluating `S_u x S_v` for a revolution. For query cylindrical coordinates `(rho, z)`, minimize the Euclidean distance to
that curve; angular minimization aligns the radial directions while the offset radius stays positive. Scalar
minimization over an enlarged parameter interval agrees with direct 3D projection to about 2e-15 on the first rejected
junction: a cross-check of that measurement, no whole-domain existence or global-minimum certificate.
Exact versus approximated offsets: [PMC09 §11.3.7]. A continuous contact construction still needs branch ownership and
endpoint compatibility; the initial-component assumptions of [Mukundan04 §1, 2.2, 4] are not supplied by an independently
continued approximate contact, and neither the reference observation nor the finite-point bound closes that problem.

### Contact ownership on a joined blend sheet

A final construction, attempted only after the previous bounded and extended
pipelines refuse, refines the **final joined sheet's** ownership intervals.
An earlier rail split is insufficient when fitting and joining change the
contact chart. Isolated original-contact/beyond-surface intersections are
retained only when they project onto the actual shared boundary edge within
the existing tolerance. They partition increasing active rail intervals;
original-face IN/ON intervals retain that owner, other intervals retain the
existing beyond-face owner. Adjacent pieces share a witnessed Split boundary.
The 3D sheet is unchanged. Copied 2D fitted contact curves are restricted by
`Segment` at interior events, retaining the original outer continuation tails.
A bounded notch crossing one far edge twice retains the two outside edge
intervals together rather than applying incompatible sequential endpoint cuts.

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

This is a numerical construction for finitely many isolated boundary events
on a regular contact strip, not a certified enumeration of every possible
contact component. Overlap segments and coincident events on distinct edges
are refused, and decreasing active charts are not refined. Midpoint ownership
classification presupposes that the witnessed boundary events delimit the
relevant owner changes; the unchanged final topology, tolerance, continuity
and collision checks remain mandatory. Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
section 5.8.2.3, equations 5.100--5.105, explain why marching requires starting
points for every component and why tangencies require additional treatment
([chapter](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html)).
The OCCT `GeomAPI_IntCS` interface separately exposes isolated points and
tangential segments; a zero point count alone is not a complete empty-set
oracle ([reference](https://occt3d.com/dev/doc/refman/html/class_geom_a_p_i___int_c_s.html)).

### Native test host reference lifetime

The native test host keeps shape entries in `std::deque`. Its `get()` interface returns a
reference; an operation can retain that reference while `store()` appends an intermediate
shape. A vector can invalidate it when it reallocates, making an otherwise identical operation
depend on the number of earlier stores. This is a host-lifetime correction, not a change to
the geometric construction. ISO C++ working draft, *Containers library*,
[deque.modifiers], paragraph 1 (accessed 2026), guarantees that end insertion preserves
references to existing elements: <https://eel.is/c++draft/deque.modifiers#1>.

### Centre-path arc coordinates for a defective sheet

A rail-parameter knot can be a tangential acceleration break rather than a
geometric curvature break. For a regular centre curve c(t), put v=|c_t| and
sigma(t)=integral v dt. The chain rule gives
`p_sigma = p_t/v` and
`p_sigma_sigma = p_tt/v^2 - p_t (c_t.c_tt)/v^4`
for each contact foot p(t). This cancellation is **derived** from the arc-length
relations in Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer
Aided Design and Manufacturing* (2009), section 2.1, equations (2.6)-(2.11),
[arc length and tangent vector](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node22.html).
It removes only tangential parameter acceleration; it cannot remove a geometric
support break, a section-policy break, a singular centre path or a varying
radius law's break. Peters, *Geometric Continuity* (2001), sections 2-3,
[reparameterisation conditions](https://www.cise.ufl.edu/research/SurfLab/papers/01handbook.pdf),
likewise distinguishes geometric regularity from equality of parameter jets.

The experimental reconstruction keeps the existing result unless a sampled
one-sided second-form comparison exceeds the existing internal curvature bound
at an inherited rail knot. For constant-size curvature-continuous sections,
centres are interpolated in the original edge coordinate, retaining its knots,
and refined to the station-solve length resolution. Numerical arc integration
and inversion then define the fit coordinate; the radius law and all station
solves continue to use the original edge coordinate. The experimental refit now
keeps a checked reconstruction even if it needs more longitudinal coefficients.
Changing the fit coordinate does not preserve the complexity of its interpolant;
coefficient count alone is not an error estimate. Position and contact checks in
`fit` remain mandatory. Additional coefficients are an explicit representation
cost, not evidence of a smaller geometric error or Pareto acceptance.
The shared rail carries
both coordinate directions, so contact, section and attachment readers agree.
The sheet fit drops only edge-coordinate knots and retains support and policy
events. Support-split rails retain their previous construction: coupling their
terminal jets requires a common group coordinate. Variable-size laws and
edge-supported contacts also retain their previous construction.

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

For a regular station solve, the implicit-function theorem transports the
regularity of its input data to the centre. A support of class C^c has normals
of class C^(c-1), so the centre carrier retains a quintic knot of multiplicity
`5-max(0,c-1)` at either path's crossing. A tangent continuation of a finite
support is only C1 in general and thus requires a C0 carrier knot. These
regularity statements are **derived**, conditional on a nonsingular station
Jacobian. Knot multiplicity and continuity follow the B-spline definition
(`Geom_BSplineCurve`, Detailed Description). The existing event locator is
shared by the section and centre fits; the latter uses first-order support
data rather than the section law's second-order data. Connected pieces map
their events into the delivering carrier's common source coordinate.

During connected centre-carrier interpolation, station-plane origins outside a parameter
curve's domain use the existing adaptor's cached polynomial evaluation. The
interpolator resolves positions to the station-solve tolerance, substantially
finer than the section fit. Ordinary station construction retains its original
parameter-curve evaluation. Each connected carrier starts from the delivering pieces'
station seeds. Single-rail centre fits retain their original event inventory and station cache.
Only a station interval extending beyond a source parameter curve's domain
needs a separate cache evaluated with the endpoint polynomial. The adaptor
evaluates the same endpoint span,
while avoiding repeated extrapolation through signed de Boor coefficients.
Horner evaluation's coefficient backward-error analysis is described by
Higham, *Accuracy and Stability of Numerical Algorithms*, second edition
(2002), section 5.1,
[Horner's method](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch5).
It is not a universal forward-accuracy guarantee for extrapolation. Independent
long-double de Boor comparisons distinguish evaluation error from a geometric
break. Neither the centre-fit length resolution nor any final geometric bound
is changed. Position sampling still does not certify the carrier's length or
derivatives; final contact, section, continuity and topology checks remain
necessary.

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

The G2 flat-support cutoff now activates a whole section profile with
`H(t)=6t^5-15t^4+10t^3`, between one and two fit-tolerance units of
`abs(curvature)*radius^2`. Both endpoint derivatives through order two vanish.
Whole-profile interpolation is necessary when a host's flat and curved laws
have different zero-curvature limits. Contact points and tangent directions are
shared by the two profiles; intermediate endpoint curvature is eased too.
The band width is a policy choice, not a new accuracy allowance. Output contact,
section, validity and collision checks remain independent requirements.
Source: Ken Perlin, *Implementing Improved Perlin Noise*, GPU Gems (2004),
section 5.4, [quintic Hermite interpolation](https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-5-implementing-improved-perlin-noise).
This source establishes the scalar interpolation regularity, not feasibility of
arbitrary G2 blend data. A host with C2 selections should return a constant
section regime: changing labels would otherwise request unnecessary C0 fit knots.

### Experimental terminal restriction continuation

The last recovery construction can continue a displaced support's terminal
contact along the identified two-offset branch. A local tangent carrier predicts
the next centre; the existing exact-offset Newton equations correct it. A sign
change against the stopping support is bisected and retained only when the exact
contact also projects onto its finite boundary arc within the unchanged bounds.
The endpoint records its support, arc, arc parameter and station parameter.
The unchanged attachment and final solid checks still decide delivery.

Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design
and Manufacturing* (2009), section 5.8.2.3, equations (5.101)-(5.105), describe
regular transverse-intersection marching and its component-seed requirement
([chapter](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html)).
Isabelle Grignon, *ChFiDS_CommonPoint* (1993), `SetArc` and `SetParameter`,
records separate restriction and spine parameters in OCCT
([source](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_CommonPoint.hxx)).
The local predictor and signed-distance boundary search here are derived;
the step size and finite exploration count are numerical policies, not existence
bounds. This WIP does not certify unsampled regularity, all boundary roots,
or consistent terminal face closure. The tested displaced contacts reach their
restriction arcs; at this stage complete solid construction still requires the
boundary arrangement described below.

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

This graph chart and knot arrangement are **derived**. They require a regular
branch that is a graph in the chosen coordinate; Newton convergence, the finite
refinement policy and midpoint checks are numerical tests, not a global existence
or approximation theorem. Tangencies, other components and branches leaving that
chart remain outside this construction. Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
sections 5.8.1.1-5.8.1.2, equations (5.84)-(5.87), reduce a patch/plane
intersection to a scalar parameter-space curve and describe coordinate stepping
with Newton correction, including the need to change charts near vertical branches
([formulation](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node101.html),
[tracing](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node102.html)).
The generated two-radius, two-scale, two-pose single/chain test checks acceptance
and the rolling-ball defect against the input skin independently of the solver.

For an offset-intersection station route, centre-length fitting retains that
route's spline carrier, including its circuit seam and orientation. It does not
re-march the source edge. The carrier is segmented at circuit wraps and affinely
mapped to source coordinates before spanwise length integration. This coordinate
rule applies to every such G2 route, independently of a near-limit knot verdict.

### Carried-sheet corner boundaries

The affine three-sheet tensor is retained when all three solved contact
crossings lie within the existing positional allowance of the incident tangent
planes spanned by the outgoing edge directions. This tests the construction
premise from geometry rather than the representation of the support. The
coordinates, poles, weights and boundary construction of that tensor are
unchanged. When a crossing is off its incident plane, the noncircular
three-sheet construction uses the intersections of the
carried contacts as its vertices, and connects them by curves in each carried
sheet's parameter domain. These curves are approximated spatially and used as
supported boundaries of the existing constrained fill. The carried sheets are
not changed. The two positional approximations share the existing allowance:
if the endpoint discrepancy is `g` and the total allowance is `B`, each receives
`(B-g)/2`. The triangle inequality then accounts for both approximations and
the endpoint mismatch; this does not certify unsampled approximation errors.

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

### Split corners of more than three conic blends

Where more than three conic blends meet at a vertex and every face there is a plane, each blend's
sheet is its section's extrusion along its edge, `S(s, r) = C(s) + r d`, with `C` one rational
quadratic. The corner is split about a centre `C` into one four-sided face at each place `Q_i` where
the contacts of blends `i` and `i + 1` cross, bounded by a half of either blend's end and by two
inner curves from the ends' splits `M_i` to `C`: the central split of an n-sided region [Charrot
and Gregory 1984; Chiyokura and Kimura 1983]. Each face is a rational Bézier patch of degree 6,
made directly; nothing is fitted.

*Ends (derived).* Blend `i`'s end runs on its sheet at `r(s) = R(s)/W(s)`, `W` the section's
denominator and `R` the cubic with the crossings' rail places and slopes such that the end leaves
`Q_(i-1)` along `d_(i-1)` and arrives at `Q_i` along `d_(i+1)`. It is a rational cubic lying on the
sheet exactly, with an exact parameter curve there. Both directions and both sheets' tangent planes
at `Q_i` are those of the face the two blends share.

*Across an end (derived).* If a patch's second row of homogeneous poles is its first moved along
`d` at the same weights, `H_j1 = H_j0 + w_j0 l_j (d, 0)`, then its derivative across that side is
`6 (sum_j B_j w_j0 l_j) / W d`, along the sheet's extrusion, so the patch is G1 to the sheet for any
offsets `l_j`. At `Q_i` the two rules give the same pole, because either end's tangent there is the
other side's direction. G2 to both sheets is impossible at `Q_i`: along the shared face the two
sheets bend differently there (the corner obstruction of [G2 at a corner](#g2-at-a-corner)).

*Halves.* Each half end is brought to the form whose two first and two last weights are 1: a
Möbius change of parameter (standard form), then a common cubic factor `(1, 2 - w_1, 2 - w_2, 1)` of
numerator and denominator, which gives a sextic. The end is split where the halves' speeds agree, so
the inner curves have unit weights and the halves join C1 at `M_i`.

*Across an inner curve.* Let `T` be the curve's tangent and `D_a`, `D_b` the derivatives across it
of the patch it ends at `u = 1` and of the one it ends at `v = 1`. The condition
`D_a + D_b = 2 cos(2 pi / n) t T` is the G1 condition of a split of `n` patches about one point
[Loop 1994; Peters 2000]. It is linear in the poles. Its first two coefficients hold by
construction: the offsets next to `M_i` are `l_i (1 - cos(2 pi / n) / 6)`, with `l_i` at `M_i`
(derived). The remaining coefficients are solved for together with the free poles, which are kept
fair by the least thin-plate energy of the nets. That is one symmetric system for all three
coordinates, with one more row that holds `C` at a point `c`; it is solved by singular values,
because the centre's equations are dependent. Its solution is `x0 + g c`, with one `g` for every
coordinate. Its last coefficient makes the inner curves' tangents at `C` an affine image of the
regular n-gon's.

*Centre.* `c` is where the surface bends least: the least bending energy `integral (k1^2 + k2^2) dA`
over the patches [Moreton and Séquin 1992]. The thin-plate energy of the nets stands in for it only
where a patch is parametrised evenly. On uneven corners the patches differ widely in shape (a narrow
face's crossing lies far from the corner), and the nets' energy pulls the centre towards the long
patches, below the corner's top, where the patch at a near crossing has to bend back sharply
(measured: at the optimum the corner bends 4 times as tightly as its blends on a member where the
nets' centre gives 130). Only `c` is chosen so; every other pole stays the nets' fairest for it, so
the conditions above hold for every `c`. The energy is a sum of squares of the shape operator's
entries in an orthonormal frame, `k1^2 + k2^2 = a^2 + 2 b^2 + c^2`, at the points of a composite
Gauss rule (4 x 4 squares of 8 x 8 points per patch: the minimiser moves by less than 0.7 % of the
corner's size from that of 8 x 8 squares of 6 x 6 points on the uneven family). It is minimised from
the nets' fairest centre, which is closed form along `x0 + g c`, by quasi-Newton steps (BFGS from
the Gauss-Newton matrix [Dennis and Schnabel 1983, 9.2]), each halved until it lowers the energy,
until none lowers it by more than its rounding. The result therefore never bends more than the
nets' centre does, by this energy.

*Orientation.* At a crossing `Q_i` a patch is tangent to the face blends `i` and `i + 1` share, and
its sides leave along both blends' edges in that face; so its normal there is that face's, up to
sign, exactly. Each patch faces out of the body as that face does. (Its normal at its middle, against
the sum of the faces' normals, decided this before; on steep uneven corners that normal is nearly
perpendicular to the sum, and a sign error refused 24 valid corners as folded.)
`l_i` is a sixth of the chord from `M_i` to the splits' mean: the first leg of a straight curve to
the centre's likely place. This is a judgement; in measurements the peak curvature of the corner
depends strongly on `l_i`. Choosing it by the same fairness instead left a split with no positive
offset on uneven pentagons.

*Measured.* Every seam is within rounding: gaps up to 1e-14 of the size and angles up to 1e-7 degrees
on regular apexes and pockets of four to six faces (rho 0.3 and 0.7). On regular apexes the corner
bends at most 1.9 times as tightly as its blends, and its bending energy is up to 7.6 % below the
nets' centre's. Uneven apexes (apex moved up to half the base radius, faces from 12 to 114 degrees at
the apex, heights 6 and 14, sizes 1e-3 to 1e3) build in all 144 members, without folds: the smallest
ratio of `|S_u x S_v|` to its largest on a patch is 0.020. The corner bends more than 10 times as
tightly as its blends on 22 of them, all on one pentagon whose faces span 12 to 38 degrees at the
apex at height 14 (at most 34 times); on the others at most 7.6 times. That pentagon's crossing on its
wide face lies near the corner's top; neither the centre alone nor centre and first legs together
(measured by direct search) bring it below 10. A template of tangents at `C` from the splits' actual
directions, with the matching non-symmetric condition `t_(i+1) + alpha_i t_(i-1) = beta_i t_i`
(`alpha_i` constant along a curve, so the ends are split at speed ratio `alpha_i`), was measured
earlier and delivered fewer corners. Choosing the centre by rules from the corner's points (along the
line from the boundary's mean to the vertex, or the point nearest the lines along which the inner
curves start), or the nets' energy weighted by the patches' side lengths or taken relative to their
Coons nets, did not lower the peak curvature on every member.

Charrot and Gregory, *A pentagonal surface patch for computer aided geometric design*, CAGD 1 (1984).
Chiyokura and Kimura, *Design of solids with free-form surfaces*, SIGGRAPH (1983).
Loop, *Smooth spline surfaces over irregular meshes*, SIGGRAPH (1994).
Peters, *Patching Catmull-Clark meshes*, SIGGRAPH (2000).
Moreton and Séquin, *Functional optimization for fair surface design*, SIGGRAPH (1992).
Dennis and Schnabel, *Numerical Methods for Unconstrained Optimization and Nonlinear Equations*,
Prentice-Hall (1983).
