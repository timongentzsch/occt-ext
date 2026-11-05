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

*Beyond the reference.* Three constructions go past what the reference modeller makes, and stay: a chamfer whose
distances vary along the edge (the reference modeller keeps the base distance or refuses), G2 where a narrow side is
lower than the ball (h < r; the reference modeller refuses it), and a chamfer larger than a round it runs over
(the reference modeller does not make it). Each is held to the same definitions and checks as the shared cases.

*Unproven.* Every marching limit (step counts, stall windows, Newton tolerances), the carry's stray bound
(4 r, or 4 / sqrt 2 times the ball's distance at the end in a sharper wedge), the radius/curvature
refusal margins, the choice between arcs and spline sections at corners, and the number of fit
refinements. There is no published bound for a G2 conic or quintic section as a function of the support's
curvature [searched; none found]; curvature-continuous requests therefore refuse in this tree until a caller
supplies a policy. Where a result is delivered it is measured again from outside (seam gaps, tangent angle,
curvature difference) and refused if it misses; that measurement is sampled, not a global certificate.

*Stock.* `BRepFilletAPI_MakeFillet` (ChFi3d). See the table below.

## BRepPatchExt: boundary-first fill

*Problem.* Given a closed loop of edges, each optionally asked to meet a neighbouring face to G0, G1 or
G2, and optional interior curves and points, build one face. The construction tries, in order, a plane,
an exact rational continuation of a torus or ball, a four-sided surface with the sides as its edges
(a Coons patch of the four curves with its inner poles moved to meet G1 and G2), and for any loop a
fitted surface over the loop's mean plane trimmed by the loop, with a net refined where the rim turns fast.
The fit minimises a weighted energy of squared second and third parameter derivatives subject to the
constraints; the normal equations are positive semi-definite for a connected net, and are solved by a skyline
Cholesky factorisation with a multiplier outer iteration. The chart over the loop uses cotangent-weight
(Dirichlet-energy) coordinates that are injective for a convex boundary [Floater 2003, Thm 4.1 (statement
read, proof not)]; for a non-convex boundary injectivity is not guaranteed [Karni, Gotsman, Gortler;
Weber and Zorin], which is why the result is looked at inside its loop before it is returned. G1 and G2
corner compatibility at a vertex where several patches meet follows the vertex-enclosure condition
[Hermann, Peters, Strotman 2009, all five pages read].

*Complete curvature normalization (derived).* A seam compares the two second
fundamental forms in one orthonormal tangent frame. For a symmetric form
`[[a,c],[c,b]]`, its spectral norm is `|a/2+b/2| + hypot(a/2-b/2,c)`.
The absolute mismatch is the norm of the difference; the relative mismatch is
that norm divided by the larger norm of the two complete forms (capped at one).
Only two exactly zero forms have zero normalization. In particular, the saddle
`z=k*x*y` at the origin has norm `|k|` even though both diagonal entries vanish.
This follows by solving the two-dimensional symmetric eigenproblem, using the
normal-curvature characterization in Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
[sections 3.3--3.4, equations (3.26)--(3.43)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node30.html).
The normalization formula is derived here; it is not an approximation-error bound.

*Parameter units (derived).* Surface-jet rank is tested on independently normalized
tangent columns. Their cross product is the sine of their angle; zero length
columns or a sine at floating-point resolution are unresolved. Coordinates use
the cross-product dual basis when the Gram determinant is poorly scaled; the
well-scaled path retains its previous arithmetic. This removes absolute UV-unit
cutoffs without asserting regularity of near-dependent columns. The numerical
resolution is eight machine epsilons, a roundoff policy, not a geometric tolerance.
The dual identities are derived from `x.(y cross (x cross y)) = |x cross y|^2`
for unit columns. The second fundamental form then acts on physical tangent
vectors as in Patrikalakis, Maekawa and Cho (2009), section 3.3 above.

*Sampled normal rotation.* Opposite normals at two separated observations trigger
midpoint subdivision rather than a fold verdict. A regular cylinder can turn by
more than a right angle between samples under a nonuniform chart. Subdivision
continues until consecutive normals agree in sign, a sampled differential loses
rank, or floating-point UV resolution is exhausted. Unresolved rank or rotation
is reported as unresolved regularity. This is a derived diagnostic, not a
certificate of local regularity everywhere or of global injectivity. The separate
Bernstein regularity certificate on chart fills remains required. The test family
uses the exact rational circle parameterization documented by Open CASCADE,
[GeomConvert::CurveToBSplineCurve, 8.0.1](https://occt3d.com/dev/doc/refman/html/class_geom_convert.html).
Automatic candidate selection retains its established preference for a graph
when a four-sided surface has coarse opposed normals. This is a preference,
not a fold verdict: the rotating candidate remains available if the preferred
construction cannot deliver, and is then checked with adaptive resolution.
This preserves previously delivered graph geometry. No interior-curvature
allowance is changed.

*Mean-plane degeneracy (derived).* Polygon area is no longer compared with a
fixed share of squared perimeter: that would impose an aspect-ratio limit on
exact planar rectangles. The centered area sum is unresolved only below an
arithmetic-scale allowance `8*n*epsilon*perimeter^2`, where `n` is its number of
vertices. Each centered radius is bounded by the perimeter, so bounding cross
products and accumulated additions gives this conservative first-order roundoff
policy. It is not a lower bound on geometric area or on the existence of a fill;
all later geometric and tolerance checks remain required.

*Guides with positional tolerance (derived).* A guide is a positional constraint,
not an exact endpoint two-jet. An endpoint tangent or normal-curvature mismatch
therefore cannot veto construction. For example `(x,y,z)=(.01*t,50,h*t)` on
`0<=t<=1` lies within `h` of the planar square for every `h<=tolerance`, regardless
of its endpoint angle. Its exact planar fill meets the positional and supported
G1 request. The exact-interpolation endpoint veto is removed; construction and
final measurements decide feasibility at the requested tolerance. This does not
assert existence for arbitrary guides.

*Interior-cell witnesses (derived).* A parameter rectangle is used as an owned
surface witness only if its center is inside the delivered trim and no boundary
arc enclosure can intersect it. Positive-weight Bezier hull tubes enclose the
owned pcurves; unresolved conversions or offset approximations yield no interior
witness. This uses the local convex-hull property in Patrikalakis, Maekawa and
Cho, *Shape Interrogation for Computer Aided Design and Manufacturing* (2009),
[section 1.4.2](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
The rectangle selection argument is derived here. It neither changes the cells
selected for surface regularity checks nor proves injectivity of a fitted face.

*Interior trim checks (derived).* Requested points are measured against the
delivered face, including its boundary. Positive-weight guide splines are
subdivided until their control hull has a tolerance-ball cover about a measured
face witness. If `d` is the witness distance and `r` bounds the hull radius,
`d+r` bounds every point by the triangle inequality. For a convex planar trim,
its tolerance neighborhood is convex, so checking every hull pole suffices.
This uses the convex-hull property in Patrikalakis, Maekawa and Cho (2009),
[section 1.4.2](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
For polynomial guide and surface pieces, an alternative paired cover avoids
subdividing to tolerance-sized balls: restrict surface cells to an in-trim
UV rectangle and follow its diagonal. Partition the certificate at the guide's
knots and at every crossing of a surface knot by the UV chord. Exact knot
insertion and restriction give polynomial intervals on which the comparison
below applies (Patrikalakis, Maekawa and Cho, 2009,
[section 1.4.3, equations 1.74-1.76](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node18.html)).
This does not force dyadic geometric subdivision to converge on
non-dyadic knot parameters. Already validated point witnesses are reused only
for exactly equal spatial coordinates on the immutable delivered face.
For bidegree `(p,q)`, its diagonal poles
are `R_k = sum_(i+j=k) C(p,i) C(q,j) P_ij / C(p+q,k)`. This follows by
multiplying Bernstein basis polynomials in the tensor-product representation
(Patrikalakis, Maekawa and Cho, 2009,
[section 1.3.6, equation 1.56](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node14.html)).
After degree elevation, the largest distance between paired guide and trace
poles bounds the whole difference curve by the convex-hull property. This
bound and the rectangle argument are derived in exact arithmetic; floating-point
conversion and pole calculations are not interval certificates. A local
projection is used only to supply an actual in-trim upper-distance witness;
failure to find one falls back to the global trimmed-face distance query.
For non-spline curves requiring approximation, the converter's reported maximum
error is reserved in the budget; that reported error is not an independently
validated interval bound. Unresolved coverage is refused, never accepted from
finite guide samples alone. Coverage can be costly on curved faces at tight
tolerances; this is an optimization item rather than a guide pole-count limit.

*Energy assembly (derived).* A quadrature cell touches local poles `(a, b)` in order
`p = (firstU + a) countV + firstV + b`. These indices increase strictly, so an element's
lower triangle is already in the global lower triangle. The assembly reads the given
pole's coordinate, or its dot product with the height direction, once per fit, and its
local-to-global indices once per element. The given net is unchanged until the solve ends.
The energy does not couple spatial coordinates, so visiting one coordinate at a time
preserves each right-hand-side entry's sequence of subtractions. Each matrix entry
likewise receives the same additions in the same order, including the original omission
of zero entries from the skyline. This removes indexing and repeated projections rather
than changing quadrature, arithmetic reductions, constraints, or refinement decisions.
An independent generated family recovers the zero-energy affine plane from perturbed
poles at degrees three and five, three scales, with separate, coupled and directional
unknowns. Four fixed corners determine the affine nullspace uniquely.

An open boundary asked for Quad sided whose runs and gaps make four sides with corners is filled by the
four-sided surface untrimmed: the poles of a gap side are free like the inner ones, and nothing is asked
along it (a side with more poles at degree 5 than a fitted net may have is first brought down to degree 5
within a tenth of the tolerance, as an arc is). Where that quad is not a face to keep, the boundary is
filled as any open one is, as before. Where the two rims face each other across two gaps, the fill must turn
from the tangent plane of the one to that of the other, often far more tightly than either rim bends (an
elbow from the end of a tube onto the end of a plate). The inside is then judged against the *bridge*
between the rims as well: along each crossing, the Hermite curve that leaves and meets each rim square to it
in its tangent plane, at the speed of the chord, with the face's bend that way where G2 is asked. It is
built from the rims' places, tangent planes and bends only, so it is pose- and scale-covariant and a bulge,
ripple or spike that leaves the rims unchanged leaves it unchanged; the fill may bend four times as tightly
as the rim or as the bridge within a tenth of the crossing's length of the same place. With no gap, or gaps
that are not opposite, the rule is the previous one.

*A section as reference.* After a guided height fit fails, a planar section guide between two
disconnected rims supplies a reference candidate. Its plane intersects each rim at an anchor.
Quintic Hermite segments extend the unchanged guide to these anchors. At a smooth rim the
extension's first and second derivatives come from the support jet, differentiated subject to
the section-plane constraint. An unsupported rim asks for position alone. The reference is the
single tensor spline `S(u,v) = R(u) + C(v) - R(u*)`, made by adding the two curves' poles and
subtracting the anchor, without a pole-count refusal.

*Derivation and scope.* This translation member of curve-network blending contains the section
algebraically. If the support is `R(u) + Q(v)` with its transverse two-jet independent of `u`,
matching the section's two-jet to `Q` at the anchor gives the same surface two-jet along the
whole supported rim. A rail whose projection derivative on the section-plane normal stays
away from zero and a simple regular planar section give a regular injective translation surface. These are sufficient
conditions for this class, not an existence theorem for arbitrary networks. Independently
parameterised opposite rims are still measured against their original curves. The candidate
must pass the original delivered-rim, support, guide, validity, interior and sewing checks;
a failed reference retains the earlier refusal. Existing deliveries are retained.

The curve-network Boolean sum and its conversion to one tensor spline are described by
Siggel, Kleinert, Stollenwerk and Maierl, *TiGL -- An Open Source Computational Geometry Library
for Parametric Aircraft Design* (2018), [sections 3.2.1--3.2.2, equations (10), (14)](https://arxiv.org/html/1810.10795v1).
Support-derived second-order ribbons and their continuity conditions are described by
Salvi and Varady, *Multi-sided Surfaces with Curvature Continuity* (2014),
[section 3.3.2 and Appendix](https://salvi.chaosnet.org/texts/ribbons.pdf).
The section-plane acceleration correction and this restricted reference route are derived here.
The interior curvature allowance remains the same unproven judgement; a compatible guide can
therefore still be refused by it.

*Planar collars and graph directions (open construction).* A planar support with unit normal
`n` is locally a height over a plane with unit normal `d` exactly when `n.d != 0`.
**Derived:** orthogonal projection restricted to its tangent plane has singular values
`1` and `|n.d|`, so its inverse condition number is `1/|n.d|`; the height gradient has norm
`sqrt(1/(n.d)^2 - 1)`. The existing `0.8` graph-selection threshold is therefore a conditioning
policy, not an existence condition. Three mutually perpendicular support normals have best
common margin `1/sqrt(3)`, which is nonzero. A simple projected boundary and a global
single-valued interior must also be established; a local inverse alone does not give them.

For a polynomial tensor graph, affine reproduction gives exact planar boundary jets by setting
every height coefficient whose closed basis support meets a planar rim to that plane's affine
height at the tensor Greville abscissae. Local support makes the value and its first two
derivatives planar at the rim; degree five and simple knots give the needed C2 continuity.
Incompatible planes must not prescribe the same coefficient. Refinement separates disjoint
collars; it cannot resolve inconsistent jets at a shared corner. The identities follow from
Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing*
(2009), [section 1.4.1, equations (1.58)--(1.61)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node16.html),
by applying polynomial reproduction in both tensor directions. The coefficient prescription
and projection singular-value calculation above are derived here.
For completeness, `xi_i = (t_(i+1)+...+t_(i+p))/p` gives
`p (xi_i-xi_(i-1))/(t_(i+p)-t_i) = 1`. Substitution in the derivative identity (1.61)
therefore gives `(sum xi_i N_i,p)' = sum N_i,p-1 = 1`; its value at the clamped left
endpoint fixes the constant, so `sum xi_i N_i,p = u`. Partition of unity then reproduces
any affine height, and tensor products reproduce its two-variable form.

This construction is not enabled. In a three-direction planar-rim experiment, 4,761, 17,689
and 68,121 coefficients met the surface-side rim checks before sewing, but the interior-curvature
ratios were approximately 15, 11 and 10 against the existing allowance of 4. These are failed
candidates, not a proof of infeasibility. Exact planar jets alone do not establish interior
fairness, topology, or agreement with the intended shape. The original routes and their
acceptance limits remain in force.

Where two sides asked for G2 meet at a corner whose faces bend differently (an arm of a round meeting a
strip whose edge bends across the arm's axis), no surface that is twice differentiable at the corner has
both bends there. It has the one bend that the two edges and the tangent planes along them leave it: four
numbers (each edge's normal curvature and the turn of its tangent plane along it) for the three of one
bend, by least squares, with each face's bend across its edge counting a thousandth as much, so that it
decides only what the four leave open (an edge that runs straight on through the corner, where the bend
across becomes the mean of the faces'). The bend each side asks across itself then eases from that one at
the corner to its face's own, as 1 - (10 t^3 - 15 t^4 + 6 t^5) of the difference over the length of the
ease: the quintic whose first two derivatives vanish at both ends, so the bend asked stays twice
differentiable along the side and agrees with the other side's to second order at the corner. Within the
ease a miss is a share of no less than the ease itself (near a corner of straight edges both the bend asked
and the surface's are near flat). The fill is held to the eased bends; the per-edge curvature it reports
is against the faces themselves, so the corner's miss is reported and the fill is not refused. Over how
long (derived here): with the faces' bends met on both sides beyond the ease, the fill's bend has to turn
across the corner by the whole difference at every distance d from it (on the one edge it is that edge's
own curvature, on the other the other face's), so it changes by that over an arc of about d; a net
follows it only where d is no shorter than its span. The ease is one span of the coarsest net (the loop's
size over 6), and at most half its side. Where no fill that eases so meets the faces beyond (three rounds
round a box corner at a tolerance of a millionth; unequal rounds and their planes round a spherical cut),
the refusal names the corner and says that G1 can be had, as before. The G2 analogue of the
vertex-enclosure condition is stated open in [Hermann, Peters, Strotman 2009]; this is the elementary part
of it that a fill can be held to.

*Corner compatibility (derived).* Take a common oriented tangent plane at a corner p.
For a regular C2 fill F, its second fundamental form II_F is continuous there. G2 contact to
support A along one incident edge gives II_F = II_A on that edge; taking its limit at p gives
II_F(p) = II_A(p). The other edge similarly requires II_F(p) = II_B(p). Thus
II_A(p) = II_B(p) is necessary. A nonzero difference of these forms rules out exact G2 at
that corner, however many poles or internally G2 patches the fill has. This is the elementary
corner obstruction, not a sufficient existence theorem for compatible data. The definition
of II and its relation to normal curvature are in Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing* (2009 online edition),
[section 3.3, equations (3.26)-(3.30), Theorem 3.3.1](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node29.html).

The scalar twist compatibility residual alone cannot detect this obstruction: for example,
in a common orthonormal frame II_A = diag(k, 0) and II_B = 0 have zero mixed term but differ
by |k| in normal curvature. This covers a bending strip meeting a flat strip, and also unequal
curvature forms of several rounds at a corner. These are impossible exact corner requests,
not missing solver capabilities. An eased fill reports its actual miss against the supports;
its ease is a relaxation of the request. Failure to meet the unchanged limits *beyond* that
ease remains a finite approximation issue, not a consequence of this nonexistence proof.
A seam with lower continuity or a singular vertex would change the regular G2 request.

*Over a chart of its own.* A fill over the mean plane is a height over it, and a height has a slope of
tan(acos |n . z|) where the rim's faces have normal n; the plane's check refuses a loop where some
|n . z| < 0.17 or where the loop seen along z crosses itself. Some loops are a height over no plane at
all: where the normals of the faces along the rim sweep a closed half circle or more (the nose of a
wing, a window round more than half a cylinder), every direction d is square to one of them by the
intermediate value theorem, so every height over every plane is vertical somewhere on its rim. Which
plane comes nearest is read exactly: with the signs s of n . d fixed near the best direction of a
search over the half sphere (a Fibonacci lattice of 4096 directions laid out in the loop's frame, so
that it turns with the loop), max over unit d of min s n . d is the distance from the origin to the
convex hull of the signed normals [Boyd and Vandenberghe, Convex Optimization, 2004, 8.1.1: the
support function], found by Gilbert's iteration [IEEE J. Robotics and Automation 2, 1966]; every
answer is a direction's own least, so a bound from below. A loop that the plane's check refuses, and
for which that margin, read at 65 places of every side asked G1 or G2 with their ends, is below the same
0.17, is filled over a chart of its own (`overChart` in `BRepPatchExt_Patch.cxx`, the chart in
`BRepPatchExt_Chart.cxx`). Every other loop is filled as before, and a loop that the chart does not
fill either is refused as it was over the plane, with the same words.
The loop is first spanned by a disc of triangles: an even lattice in the unit disc, the loop on its rim
at even steps of its length (160 places round), whose places in space minimise the sum of the squared
umbrella Laplacians at every node of the rim and inside it, the rim held (discrete thin-plate fairing
[Kobbelt, Discrete fairing, 1997; Botsch and Kobbelt, SIGGRAPH 2004]). Along a side to be tangent to a
face, a ring of nodes outside the disc stands in the face, at each piece's middle less the way the fill
leaves the side (away from the face, read in the face's own parameters) times the height of an even
triangle; with the rim's umbrella near nothing, the nodes inside next to it stand where the face's
tangent plane goes on, the discrete form of the condition on the first derivative across the boundary
that two fixed rings give the biharmonic equation in Botsch and Kobbelt. Elsewhere the outer node is
free, a natural end. The umbrella is the Laplacian of the lattice's own parameters, so the disc is laid
out again, twice, as an even lattice in the chart it lays flat to, where an even triangle is near even
on the surface and the umbrella near its cotangent Laplacian [Pinkall and Polthier 1993], and faired
again; from the first of these on it is pinned to 33 places of every guide where it was nearest them
(each pin counting 100 umbrellas, a judgement: the wing tip's disc then passes within 0.006 of its guide).
It is laid flat by Boundary First Flattening with the scale held at none along the rim [Sawhney and
Crane 2017, sections 5 and 6.1; Springborn, Schroeder and Pinkall 2008, appendix E]: angles kept
everywhere and lengths along the rim, so that the thin-plate energy of the parameters is near the
surface's bending energy, which it is not over a plane where the surface stands steep. The flattening
has no injectivity guarantee for a boundary that is not convex [Floater 2003, theorem 4.1; Sawhney and
Crane, section 7]: a triangle laid flat folded over refuses the chart. A net free in space is then fitted
over the chart as over the plane, each side at its places in the chart, its tangent plane asked as
N . S_u = N . S_v = 0 and its bend as n . S_nn = II(S_n, S_n) taken about the surface of the pass before
(both conditions on the derivatives in space, none on how the chart lies), and trimmed by the sides'
curves in its parameters. A net free in space can meet every tangent plane by not leaving the rim at
all, S_n running to nothing, which is no parametrisation: measured on the wing tip without it, the
tangent planes were missed by 20 to 80 degrees with |S_v| from 1 to 120 along one side. In a chart that
keeps the rim's lengths a surface the chart is conformal to leaves the rim as fast as it runs along it,
|S_n| = |S_t| = 1, and along each side asked G1 or G2 the speed along the way the fill leaves is held
within 0.5 of that: read, and asked only where it leaves that band, an inequality and no wish. The band
is a judgement (the chart is the disc's, not the fill's). The net starts as the coarsest of the ladder
and spans are halved where a place is missed, as for the Adjacent fit, up to the poles a net free in
space is let have (`freePoles`, 640); once a net meets the limits the next finer one is fitted and the
coarser kept where they lie within `settled` of each other, the stop above. The surface is kept only
where S_u x S_v is shown away from zero on every cell the face uses: per cell it is a polynomial of
degrees (2p - 1, 2q - 1) whose Bernstein coefficients come from those of S_u and S_v [Farin, Curves and
Surfaces for CAGD, 2002, 5.8], and by the convex hull property a cell whose coefficients all have a
positive part along one direction has no zero; a cell that does not is cut in four, at most four times.
That is the surface form of the Bernstein test of a Jacobian's sign in planar parametrisation [Gravesen,
Evgrafov, Nguyen and Nortoft 2014].
Windows round a cylinder of radius 15 between heights 8 and 20, G2 to four strips of it (the application's fixture
at 200 degrees; 180, 240 and 300 held out from the design, and 270 on a cylinder of radius 20 between 10
and 34) are filled on the same nets at rest, turned, mirrored and turned at a thousandth and a thousand
times the size; read by projection onto the strips they meet every edge, the strips' tangent planes
within 0.005 degrees and their bends within a share of 0.004 (95th percentile), and stand within 0.009 of
the cylinder (the 200-degree window 0.0024, on 11 x 8 poles, in 0.34 s of CPU time natively; 300 degrees
1.4 s).
Before fitting, `netSpacing` predicts a chart's cost by replacing higher derivatives by those
of the osculating circle: a circle of curvature k has arc-length derivative magnitude k^p at
order p+1. This substitution is exact for the circle. Curvature at a point of a general curve
or surface does not bound its higher derivatives or those of its chart. Consequently the
forecast and its 640-pole cutoff are computation policy, not an approximation-error theorem,
a proof of a necessary pole count, or a proof that a refused fill does not exist.

Sande, Manni and Speleers, *Explicit error estimates for spline approximation of arbitrary
smoothness in isogeometric analysis* (2020), [section 1.2, equations (2)--(3), sections 6--7](https://arxiv.org/html/1909.03559v3),
give upper error estimates using the actual appropriate-order Sobolev seminorm and, for mapped
geometry, derivatives of the map. Their L2/Ritz estimates are not pointwise geometric curvature
bounds and cannot be inverted into a lower bound on a feasible pole count. The production
forecast still uses the circle surrogate; replacing it requires measured constructions and
unchanged final geometry checks.

An unconstrained spatial fit currently enforces one transverse curvature equation, with five
outer linearizations. G0 and G1 tolerance do not bound the omitted tangent and mixed second
fundamental-form components. A complete correction must check the actual nonlinear form and
the delivered parameter curves, including places between its fitting samples, and establish
regularity and interior acceptance. Support-derived second-order ribbons provide local jets
(Salvi and Varady 2014, section 3.3.2 above), not a convergence theorem for a single fitted
spline. Globalized constrained iterations and feasibility restoration are described by
Nocedal and Wright, *Numerical Optimization*, second edition (2006), section 18.4,
[Algorithm 18.3](https://mcube.lab.nycu.edu.tw/~cfung/docs/books/nocedal2006.pdf); their hypotheses and actual convergence must be verified for the assembled
problem. Merely adding three frozen curvature rows does not establish either.

Open:
a constrained feasibility solve with complete curvature conditions and measured convergence,
then refinement where the verified geometric residual requires it. The cap of a cone standing 80 degrees or more against its circle, refused as too steep
before, asks for 529 to 576 poles at a tolerance of 1e-4, within the 640: it is tried, and its coarse nets
are not shown regular, so that it is refused as before at about 0.4 s more per loop (the steep
cones' GTest 23.4 s against 17.8 before there was a chart). That is not foreseen from the chart: the
chart's own map is regular (a triangle laid flat folded over refuses it), and what is not shown
regular is S_u x S_v of the net fitted over it, which the chart does not decide.
The constants of the chart that are not derived: 160 places round the loop, 33 pins per guide, two
relays, the pin's 100, the band's 0.5, four cuts of a cell, the 4096 directions of the search, which
only seed the exact step, and the 256 bins over which the spans asked are read across the chart; the
Favard constant is taken at its upper end, so that the spans are, as far as the bound goes, enough.

*Unproven.* The interior acceptance gate (a fault when the interior bends more than 4 times the rim's
bend, with a floor of 2/size) has no derivation; a bound of the form K sqrt(2D) in the energy gap D was
tried and rejected. The bridge only replaces the reference that the 4 multiplies; the 4 and the tenth are
judgements. Over 76 tube-to-plate and wall-to-plate elbows (three heights, tilts and widths, at rest and
turned) that met their rims, fills bent at most 2.1 times their reference (up to 18 times the rims'
bend, the old reference, which refuses 36 of them). On the same elbows with injected defects that leave the
rims unchanged, every ripple of a hundredth of the size (6 or 12 waves) and 80 % of bulges of a twentieth
are refused, but gentler defects where the elbow turns pass more often than against the rims alone: on the
40 elbows the old rule keeps, bulges of a twentieth are refused 108 times in 120 instead of 120, spikes of a
fiftieth 104 instead of 120. The flat G2 limit, the rules by which a ladder that has met nothing yet gives up
(`stillFalling`, `same`, `unruly`, which compare misses of unconverged nets and so depend on rounding and
therefore on pose), the span and budget constants of the net, the `settled` share below and the 0.17 slope
limit are unproven.

*Stop.* A net that meets every limit is not yet the fill asked for: that is the fairest surface of all that
meet them, which the nets of the ladder come nearer as they become finer.
Once a net meets the limits, the ladder goes on until two nets in a row lie within `settled` = 1e-3 of the
loop's size of each other, read from twice as many places as the finer net has spans each way inside the
loop to the nearest place of the other surface; the coarser of the two that meets the limits is kept, and
where the ladder or the nets allowed end first, the finest that meets them. That is a Cauchy test: where
the nets converge at least twice as fast as the ladder refines (a ratio rho <= 1/2), the kept net is within
d / (1 - rho) <= 2 d of their limit, d the distance between the two. A share of the size makes the stop
the same at every size and pose. A net free in space and a height over the plane minimise different
energies and come nearer different surfaces: they are not one sequence, and once nets free in space meet
the limits, no height is tried after them. The rules this replaces took the first net that met the limits
(closed loops) or a coarser net unless a finer one bent less by 0.75 (open loops); measured on the
kinds of published example, that delivered a seat 9.2 off its master (limit 1.04; every finer net 0.17), a mouse top
0.61 (limit 0.6; finer 0.19), a high mouse top 3.3 (limit 0.9; 0.36), a hull bottom 3.6 at rest and 0.04
turned (limit 1.44; 0.03), with the outcome changing with rotation and translation. With the stop each is within its limit at rest, turned, at five
times the size, rotated about z and moved a thousand away, with the same reading to three or four digits.
On parameter variants held out from that design (hulls, mouse tops, seats of other proportions; three
poses) the old rules missed 17 of 36 and the stop none. Fills of polygons of three to seven sides at G1
and G2 on strips of biharmonic masters (a dome, saddle, ridge, bowl and monkey saddle) come out on the same
nets and surfaces as before. `settled` is a judgement (about the height by which a tilt of 0.1 degree moves a
surface over half the loop's size): over the nets of the native tests these distances spread from 1e-5
to 1 of the size with no gap between settled and unsettled pairs, and the seat, mouse tops and hull
bottom above are delivered on the same nets for any value from 3.1e-4 to 2.4e-3. The cost is at least one net beyond the one kept: about
2.5 times the CPU time of a fill whose first net was already settled, 5 times for the seat, 14 for the high
mouse top, and nets that a coarse net met only roughly come out finer (a corner of rounds 2 / 3 / 5: 43 x
48 poles before, 82 x 89 now; the coarse net was 3.5e-3 of the size from the next). Where nets keep
differing by more than `settled` the whole ladder is run (a G0 fill at flatness 0.95: 1.2e-3 to 3.8e-3
between nets, up to the finest). Open: the Cauchy bound needs the distances to fall, and they need not.
On the saved eight-edge opening (G0, Natural, a loop of 38.3) the nets lie 2.38e-3, 7.0e-4 and 7.5e-4 of the
size apart (6 x 3 to 10 x 4, to 16 x 6, to 24 x 10): the ratio is 1.07, not at most 1/2, and the 2 d bound
does not hold. The nets stand on a plateau below `settled` (their areas 70.20 and 69.98), and on such
loops the stop is a judgement, not a bound.

*A ladder that does not settle.* The same opening as a Patch at G0 never settles: its nets lie 3.8e-3,
3.1e-3, 1.3e-3, 2.9e-3 and 1.2e-3 of the size apart, up to the finest. Why they need not come nearer: a
condition is met once it is within its reach, and the multipliers stop at the first shifts that put every
condition there ("nearer is not better"); the correction then moves the net onto the curves by the least
movement of its poles. What a net delivers is so one surface of the band of all that meet the limits, the
one at which that iteration and that correction stopped, and not the least bending of the band. The band
has an inside, the limits being inequalities, and as the nets become finer there is nothing that draws
their surfaces to one point of it. Read inside the loop, the nets of this opening bend 13.8, 15.6, 7.7,
7.8, 12.2 (missing the rim, 0.0105 for 0.01) and 14.9, ∫ (k1² + k2²) dA. The rule that kept the first
net that met the limits delivered the 13.8; keeping the finest delivered the 14.9.
What is asked is the fairest surface that meets the limits, the least of that energy over the band. Every
net that meets them is a surface of the band; where the ladder has not settled, nothing says that one of
them is nearer the least than the others but their energy, and the one that bends least is, by the energy
of which the least is asked. So where no two nets in a row have come within `settled` of each other, of
the nets that met every limit the one of least bending is kept (`bendingInside`, read at the places the
inside of a face is looked at; k² dA has no unit, so the choice is the same at every size and pose to
rounding). A ladder that settles keeps what it kept. Nothing new is tuned: the energy is the measure of
fairness, and the places are those of the inside check.
Pose dependence is measured, not removed: 47 sampled breaches over 154 rigid-motion variants, none changing
whether a case is refused. The ease of one coarse span is a judgement bounded by measurement, not an
optimum: on the teddy (arms and strips of radius 50, sides 79 long, the span 24) eases of 0.3 (the
distance within which a bend that differs by as much moves a surface by no more than the tolerance), 3.9
and 7.8 left every net short of the bends asked, or meeting them and bulging inside; 16 was met on the third
net, 24 on the first to third (1.2 to 5.3 s). G2 is then met on the middle 40 % of a teddy side and on
less of a shorter one; the application's reading of a teddy seam counts the eases as misses above its 95th
percentile. Coarse nets that meet an eased loop can be lopsided: the turned teddy is filled 5 off its
three-fold symmetry on the first net, 0.4 on the third at rest. Straight-sided loops of three to six sides on
rounds tangent to their plane (corners of 60 to 120 degrees, radii 15 to 80, rising and falling, at rest and
turned and scaled, one side split between two radii) are all met. On a family held out from the design
(irregular polygons of three to six sides lifted onto a dome, a saddle and a plane, each side held to a
strip tangent to that surface and bending 0.05, -0.03, 0, 0.02, ... across it; at rest, turned and
mirrored, and at a thousandth and a thousand times the size with the tolerance in step) 27 of 36 per pose
or scale are met, with the same cases met at every pose and scale, G1 within 0.015 degrees and the strips'
bends within 1.1 % on the middle of every edge. The 9 others are the polygons of three, five and six sides
on the plane: where a strip that does not bend meets one that does, the flat strip is held beyond its
ease to the flat G2 limit, which no eased fill meets, and they are refused naming the corner. A commercial fill's
own manual says the same of such corners ("the intersection of two curves or edges determines the
generated surface's tangent plane and curvature at that point, which may not agree to the existing
surfaces'"; "disconnecting the boundary curves is a solution for enforcing G2"), and its examples
deliver with such corners marked as conflicting.

*Read up to the corners.* A fill is measured at the 129 places of Measure Continuity along each side
and, between the first and the last of them and the side's ends, at places that halve their way to
each end, down to `cornerInset` (1e-6 of the side, a fraction and so the same at every size; a
judgement of how near is at the corner: read down to 1e-5 or to 1e-7 instead, the application's Patch
tests and the mixed records came out the same). G0 and G1 are read there as everywhere; the bend is
read there too, but not within an ease the corner's faces are shown to ask for (see agreeAtCorners),
where it is read from the first place on as the ease's allowance was set. Read so densely inside the
first 1/258 of each side, fills delivered before were not all within: a five-corner tangent chain at
1e-6 left its rim by 7.97 times the tolerance at 1/16 of the first place, and a later panel at G2 bent
1.125 times its share there; a fill of the band fit at a hundred times the size of an unresolved
corner was within.

Where two sides that are only to be met meet at a corner p, a surface that holds both is tangent at p
to the plane of their two tangents: each edge lies on it, so each tangent lies in its tangent plane,
and two tangents that are not parallel span it. The plane may stand far from the surface around it
(at that chain's corner the edges meet at 1.72 degrees in plan with slopes of -0.043 and -0.016, and
the plane of the two is 42 degrees from the mean plane), but it is asked at p only: along each edge the
surface's tangent plane need hold only that edge's tangent, and turns as fast as it will. Where the
edges run on within the G1 angle the two tangents span no plane, and nothing is asked. Within the
tolerance a surface need not be tangent to it, and asking it of every such corner bends fills that
were within before (the saved eight-edge opening came out 68.66 for 70.20): so a correction onto the
sides holds, besides its even places, those of the measure's places near such a corner that the
corrected surface misses by more than the tolerance, and finer nets of the ladder are asked the plane
at that corner. On the tangent chain the surface's plane at p came within 0.04 to 0.3 degrees of it and
the miss near the corner from 7.97 to 0.27 times the tolerance; it is filled on the same net at rest,
moved and turned, mirrored, and at a thousandth and a thousand times its size (128 x 99 poles at 1e-6,
75 x 54 at 1e-5), the posed fills within 2.5e-4 of the tolerance of the one at rest mapped back.

*Size.* Lengths are shares of the loop's size or of the tolerance, with one exception by design: where a
fill goes over into planes and other faces at once (or into a cone alone) its rim and its edges' parameter
curves are made to a tenth of `Precision::Confusion()`, OpenCascade's absolute least distance, so that its
edges keep that as their tolerance and STEP gives the joint back (a G1 cap over a circle on a cone, tolerance
1e-4 in step with the size, is then made on 11 x 11 poles at 1e-3 and on 17 x 17 at 1 and 1e3). A fitted net moved onto its curves ends
its ladder when the bounded rim correction leaves its own places further off than the tolerance (the
correction's own residual, no constant). The fills of the tests, with every length and the tolerance
scaled by 1e3, are made or refused alike with the same nets; at 1e-3 so are a disc, a cap, an open loop of
arcs, the dish at G0 and G1, the square on rounds and the elbow. The dish at G0 was refused at 1e-3
as crossing itself: two pieces of one straight side, collinear, were taken to cross by the signs of
orientation determinants whose true value is nothing and which rounding alone decided, otherwise at
one size and not at another. A point now counts as on one side of a segment only where the
determinant exceeds the bound on what rounding its products and their difference can do (Shewchuk
1997, the ccwerrboundA of orient2d), a bound relative to the numbers themselves and so the same at
every size. The square on rounds, which missed G2 by its quad at 1e-3, is made alike there now; it
changed with the stop by convergence or the corner ease in the fit's own units, and was not traced
further. The pentagon and the split side took other nets at 1e-3: at a corner where both edges run
straight on rounds that are eased to one bend, everything the corner's four numbers ask is rounding,
and the share by which they disagree, rounding over rounding, counted where it came above an absolute
1e-15, at one size and not at another. It now counts only above the rounding of the faces' own
curvatures, 2^-40 of the larger of their principal curvatures (`jetRounding`, a judgement: a few
thousand units in the last place, far below any curvature asked). Its margin, as the corner's
disagreement over the larger curvature: rounding alone, on the pentagon and the split side, came to
at most 9e-17, ten thousand times below the floor (9.1e-13); the least genuine disagreement, on the
saved opening below, 6.7e-12, only 7.4 times above it. That margin is thin: a genuine
disagreement eight times smaller than this one would be taken for rounding, and its turns dropped.
A floor derived from the faces' own representation instead (each derivative off by up to
(gamma_n + eps) of the sum of |N_i^(a)| |N_j^(b)| |P_ij| over the poles where the face stands, carried
into the normal curvature) was measured and is not used. On the held-out flat strips it is about
2.7e-13 (2^-40 gives 2.7e-14 to 4.5e-14), and the largest of a corner's four numbers is:

| pose | largest of the four | derived floor |
| --- | --- | --- |
| at rest | 2.6e-19 to 3.2e-18 | 2.6e-13 to 2.7e-13 |
| mirrored | 8.4e-12 to 1.7e-11 | 2.5e-13 to 2.8e-13 |
| moved | 1.2e-12 to 1.4e-12 | 1.9e-13 to 3.7e-13 |
| turned | 2.3e-11 to 4.9e-11 | 2.8e-13 to 6.7e-13 |

The strips are interpolated from posed points, and an interpolation is rounded by the solve's
condition, not by the poles' storage: posed, they bend by about 1e-11 along their edges and are other
surfaces than at rest, which no bound from their representation can tell. And the saved opening's
corners ask turns of numbers of 1.1e-18 to 2.2e-18 (its largest curvature 4e-8), five orders below
what its representation is known to: its fill is delivered only through turns that its input cannot
know. 2^-40 stays, a judgement. Open: the flat strips come out on other nets mirrored, moved and
turned (their inputs differ); and the saved opening's real defect is the fold of the net over the
plane without those turns. The structural path is a floor at the flat share of the G2 limit, the
measure's own, once that fold is fixed. A floor at the flat share of the
G2 limit instead, as for a miss, would be the measure's own, but a saved opening of near-flat spline
faces (curvatures of 6e-7 on a loop of 38) is filled only when the turns that their small but
genuine disagreement asks are kept: without them its fill folds over. That fold is an open defect of
the net over the plane, not of the corner: a fill whose rim asks no turns should not fold, and the
turns only hide it. Open: at 1e-3 the corner of three rounds
(its tolerance below `Precision::Confusion()`) is made; the teddy is refused at 1e3. A
G1 cap over a circle on a cone is refused as too steep from the same tilt on (between 80 and 80.4
degrees) at rest, turned and at 1e-3 and 1e3, with the same nets below it; but the switch from the height
graph to the spatial net at a support normal of 0.8 against the mean plane changes the cap's area by 2.1 %
between 36.8 and 36.95 degrees (a cone as a spline), and the analytic cone's cap changes between two shapes
from 70 to 74 degrees (area 2.45, 1.23, 1.22, 1.22, 2.56 of the disc's).

*Stock.* `BRepFill_Filling` and `GeomPlate_BuildPlateSurface`: a plate fit that accepts G0, G1 and G2 point
and curve constraints and minimises a different energy. See the table below.

## BRepCurveExt and BRepSurfaceExt

*Curves.* Offsets of planar curves with three corner treatments (arc, chord, extension to the crossing) by
working on a copy moved to the origin and scaled, so that OCCT's absolute 1e-7 lengths do not change
the result with the curve's size; the outline of a body along a direction from edges and the turning curve
n.v = 0 of curved faces, with the places where visibility can change found exactly; tangent circles by Newton
in a frame of the curves' own extent; projection, rebuild with a prescribed degree and count, refinement and
regions of coplanar curves. A spiral is the helix (angle 2 pi t, radius R + p t, height h t / turns).

*Surfaces.* Alignment of the control row beside a shared boundary on the line through the boundary
pole and the neighbour's next pole, so that the control hulls make one grid; orthogonal projection of pole
offsets onto the affine constraint space; polyspline development of an oriented quad graph in integer UV;
face rebuild by interpolation checked on an independent grid; the number of poles a deformation target needs
from the C2 cubic interpolation error 5 h^4 max|f| / 384 and the end-slope term h^3 max|f| / 4;
face-separation proofs by Bernstein bounds and, for a plane beside a surface of revolution about its normal,
in closed form (the surface meets the plane in circles; a band of width gap / sin t).

### Alignment input-shape prior

`alignSurface` accepts `inputShapeInfluence` in the closed interval [0, 1], finite,
with default zero. It acts only when extra blend rows are requested. Zero influence,
or zero extra rows, executes the previous solve without adding any matrix entries.
The moving surface and reference are copied; only a certified result is stored.

The construction is Tikhonov regularisation of the free poles around their original
positions, following Akil Narayan, *Regularization and relaxation* (University of Utah,
22 February 2022), slides 6–7, “Tikhonov regularization” and “Tikhonov regularization
for least squares” ([lecture](https://www.sci.utah.edu/~akil/docs/courses/2022spring/math6880/slides09.pdf)).
Let P0 be the pole vector after exact degree/knot refinement and before alignment.
Let F(P) be the existing bending-plus-variation quadratic, with scalar-coordinate
matrix H on the n free poles. The objective is

    F(P) + w a sum_free ||P_i - P0_i||^2,  a = trace(H) / n,

subject to the existing continuity conditions and fixed remaining poles. Trace
normalisation and the [0, 1] control domain are design choices derived here, not a
published surface-specific weighting law. The trace is computed before constraint
rows or the prior are added. This makes w dimensionless and preserves rigid-motion
and uniform geometric scale covariance for a fixed basis. It does not make the
pole metric invariant to changing basis or arbitrary UV reparameterisation.

For displacement unknowns d from the current iteration Pc, add w a to each free
coordinate's diagonal and w a (P0 - Pc) to its right-hand side. Keep P0 unchanged
through all nonlinear continuity iterations. This costs O(N) time and O(N) saved
coordinates for N total poles, adds no poles, and keeps the matrix band unchanged. The generic fit
applies this term after its optional energy cache, so no prior is cached accidentally.
The unchanged condition/multiplier solve and final seam check still apply; the
prior cannot authorise a seam miss. If a surface already meets the requested
continuity, the existing unchanged-surface path still returns it.

**Derived bounds.** H is positive semidefinite, since it integrates squared spline
derivatives. For a > 0 and w > 0, H + w a I is positive definite on the free
coordinates. Thus a nonempty affine constraint set has a unique minimiser.
For minimisers P1 and P2 at w1 < w2 on the same feasible set, adding their two
optimality inequalities gives (w2 - w1)(D(P2) - D(P1)) <= 0, where D is the
squared pole displacement. Thus total displacement is nonincreasing; individual
poles or curvature need not vary monotonically. Nonlinear G2/skew iterations have
changing linearised feasible sets: the inequality proves each fixed-row problem,
not global convergence or a universal monotonicity theorem for those iterations.
The final seam certification remains sampled, not a whole-domain Hausdorff proof.

*Unproven.* Work budgets and halving depths of the separation proofs, the crease threshold, interpolation
grid sizes, and the choice of pole counts "twenty-one for any target". The sampled checks are not Hausdorff
certificates.

*Stock.* Offsets: `BRepOffsetAPI_MakeOffset`; outline: `HLRBRep` (its documented limits on tangent and
edge-on faces are what the exact event search avoids); the rest has no direct stock counterpart.

## BRepOffsetExt, BRepSkinExt, BRepBooleanExt

*Offsets and shells.* A solid offset first classifies faces the distance would shrink to nothing, removes
them from a working copy by extending their neighbours (`BRepAlgoAPI_Defeaturing`) and offsets that copy;
a face too tightly bent in part is refused by name. Stock `BRepOffset_MakeOffset` removes vanishing faces
only where every face is planar (`myIsPlanar`); a fillet whose radius the distance reaches comes back as an
invalid solid or not at all. Push-pull follows Zou and Feng (2019): the solid swept by the moved faces is built
and joined to or cut from the body by a Boolean, which resolves whatever the move does to the rest. Thicken
divides a C0 sheet's supports without changing their geometry before offsetting.

*Face offset that keeps the topology* (`BRepOffsetExt_Kept.cxx`, tried first by push-pull where every
face moves; it holds for any set of faces with no tangent neighbour that stays, but a move of some faces
keeps the engine's route until a difference it caused in a later incremental display mesh is explained). The reference modeller's face offset keeps every face (the same face
count before and after in 536 of 537 results), which is the intersection join. Each face goes onto its
offset surface, which `BRepOffset::Surface` gives in the surface's own parameters (a plane translated; a
cylinder, cone, sphere or torus with the same axes and (u, v), since S(u, v) + d n(u, v) is again such a
surface in those parameters; otherwise `Geom_OffsetSurface`). *Tangent edges:* if two faces have the same
normal n along an edge c and move by the same d, the point c + d n lies on both offset surfaces with the
normal n on both, so they touch along the offset of the edge and its curves on the faces stay; this is
tested, not assumed (the two faces' points of each edge parameter after the move agree within the
tolerance, or within the edge's own tolerance where that is larger). Two faces tangent only to an angle
θ come apart by |d| · 2 sin(θ/2) ≤ |d| θ more when both move by d. An input that is G1 only to θ, its edge
tolerance that wide (the reference modeller's G2 blends meet a cylinder to 6e-5 rad, edge tolerances 4e-5 to 2e-4),
has its edge kept, and the result's edge carries the input's tolerance plus |d| θ_max. θ_max is found,
not bounded (judgement): 8 samples in each piece of the edge where its curves on the faces are polynomial
and cross no knot line of either surface, each local maximum refined by golden-section search. A bound
would need the faces' curvatures bounded along the edge, and the hull bound on curvature was three orders
loose on blend sheets (above). Measured against 4 097 samples on every nearly tangent edge of the 22
recorded rows, it found the larger maximum on all of them; 33 even samples fell short by up to 36 per
cent. The certificate holds that edge to its carried tolerance and its vertices to the input vertex's
plus the same spread; the result is never claimed more accurate than its input. Tested on a held-out family:
prisms whose round meets its sides at θ from 1e-5 to 1e-3, at three scales. *Sharp edges:* the new edge is the crossing of the two offset surfaces, solved in the plane
through c(t) normal to the old edge (four equations S1(u1, v1) = S2(u2, v2), (S1 − c(t))·c'(t) = 0 by
Newton from the old edge's curves on the faces), so it keeps the old parameter t; the curve and the
curves on the faces are cubic interpolants of those sections, refined locally (a midpoint added where an
interpolant is further than a tenth of the tolerance from the section there). An offset surface is one
order less smooth than its basis (its normal takes a derivative), so where the edge crosses a line along
which a moved surface is less than C2, or the old edge (whose tangent sets the section planes) is less than
C3, the curve is interpolated apart on each side and joined C0 there; the crossing is found among samples
and by bisection. Newton steps that do not reduce the residual are halved; where two surfaces only touch
and the residual stalls, a point within half the tolerance of both is accepted (judgement: the other half is left to the interpolation, a tenth, and the input's own edges). *Vertices:* on a tangent edge's (or seam's) moved curve, crossed with the face
around the vertex that meets it most steeply among those that meet another at a sharp edge there (all
tangent: the vertex moves along the common normal); between two edges of the same two faces, in the
section of one of them; else where three faces cross (the three whose normals are furthest from a common
plane). A vertex that is then not on every face around it means the move changes the topology (four
planes of an apex; a face that shrinks away), and the construction declines. *Certificate,* apart from
the construction: a valid solid; points inside every face and along every edge at the moved distance from
their sources' unbounded surfaces, on the side moved to (foot point by Newton, Gauss–Newton where the
Hessian is not positive definite); faces clear of each other (below); volume grown by an outward move and shrunk
by an inward one; and where every face moves, nothing nearer the body than the distance. A failure of any
of these falls back to the routes below. Judgement: the refinement of an edge stops at 16 halvings of the first spacing or 4096 samples (a curve a cubic cannot follow there has a kink).

*Faces clear of each other.* Always the general self-interference check (`BOPAlgo_ArgumentAnalyzer`) on
the exact faces, on what the separation pre-pass (`BRepSurfaceExt_FaceSeparation.hxx`, pole hulls and
Bernstein bounds) cannot prove clear. The pre-pass does not read offset surfaces, so it runs on proxies: a
copy whose offset faces are splines within ε = a hundredth of the tolerance of them, each face's tolerance
grown by ε. Two proxies it proves further apart than their tolerances have their exact faces apart too;
that a proxy is one to one says nothing of its exact face, so each offset face also goes to the check by
itself unless the pre-pass leaves it in a part that is checked. The result keeps its exact faces. Tried and
dropped: skipping the check by the reach [Federer 1959, Thm 4.8] (|d| κ_max < 1 per face, faces that share
no vertex further apart than both moved faces reach). A rigorous κ bound per rational Bezier patch was
sound on the matrix's bodies but three orders above the sampled curvature on blend sheets, and the reach
leaves faces that share a vertex unproven.

*Offsets that load again.* The binary BREP reader builds an offset surface with OpenCascade's continuity
check (`BinTools_SurfaceSet`), which refuses a basis that is C0 at a knot unless it finds it G1 there; a
blend's sheet whose parts are joined C0 in their parameters is refused, and a body holding such an offset
cannot be loaded from a saved file. Such a face's offset is built unchecked for the construction and, once
the face's extent is known, replaced by a spline over its box within a tenth of the tolerance
(`GeomConvert_ApproxSurface`), its edges' curves, vertices' points and continuity records moved onto it
and their tolerances grown by the error; every result loads again (checked on the matrix's bodies). An edge
whose old curve is a line or a circle keeps one: the offset of a line on a plane or cylinder is a line, of a
circle round a surface of revolution's axis a coaxial circle, with the parameter kept (checked at 17
points). At a pole, where S_u × S_v vanishes against the derivatives' size, the certificate checks the
distance alone (judgement: below a millionth of |S_u|² + |S_v|²), as the normal's sign is not defined there.

*Read back as written.* A body read from a saved file has the directions it was written with, bit for bit
(`gp_Dir::Restored` in the BinTools and GeomTools readers; corrections.md, gp_Dir-restored), so a body
and the same body loaded again mesh alike. Open: the X and Y directions of an axis placement (`gp_Ax2`,
`gp_Ax3`, `gp_Ax22d`: planes, conics, elementary surfaces) are rebuilt from the main and X directions by
cross products when read, and come back exact only where those products are; a reader that keeps
orthonormal axes as written is not yet made.

*Skins and sweeps.* Frames along a path are rotation-minimising by the double reflection method [Wang,
Jüttler, Zheng, Liu 2008], with a fourth-order error in the step; stock `BRepOffsetAPI_MakePipeShell`
twists a profile only where the path starts square to it. Gordon surfaces follow [Siggel et al., section 3.2]
with compatibility by common reparametrisation; a skin's cross-section lines are natural C2 cubics, whose
minimum of the integral of |s''|^2 over H2 interpolants is proved in the utilities section above. Closed
sections are made fit by cyclic seam alignment: the cyclic cross-correlation of two sampled rings is maximised
over start shifts and orientation [Dogan, Bernal, Hagwood 2025, rotation by Kabsch/SVD proof read; the
version here fixes rotation]. Draft of walls that `BRepOffsetAPI_DraftAngle` does not take (it handles planes,
cylinders and cones) uses the ruled isocline surface over the wall's neutral curve, whose every tangent plane
makes the draft angle with the pull direction.

Where a plate filling cannot meet the constraints, a guided skin through three or more profiles with
G1 or G2 at an end is the Gordon skin of the net plus a correction across that end. Two surfaces that share an edge and their tangent planes along it meet G2 when
their normal curvatures agree in one direction across the edge (the linkage curve theorem [Pegna, Wolter
1992, J. Mech. Design 114]); the skin already holds the profile, so G1 needs only the way across, S_v, turned into the support's
tangent plane, and G2 then only the normal part of S_vv set to the support's second fundamental form at
S_v. Both moves are along the support's normal. The correction is h1(v) m1(u) + h2(v) m2(u), with h1 and h2
the quintic Hermite functions of the end (h1' = 1, h2'' = 1 there, the rest nothing) that vanish to the second
derivative at the nearest place of the next profile, so the other profiles are kept and the skin stays C2;
m1 and m2 are the moves less their straight interpolation between the guides, so the guides are kept too.
Where a guide leaves the end profile out of the support's plane, or with another bend, the surface follows
the guide there, as between two profiles. The plate filling misses the bend limit or the boundary on 22 of 24 nets of a seeded family. Successful
plate fillings are retained: exact net interpolation does not imply a smaller integral of squared principal
curvatures. On one accepted net the correction raised that integral from 0.9364 to 1.2673. Nets whose
skin carries weights remain outside this correction.

Approximating an individual conic by non-rational cubic pieces is not the unresolved problem of
correspondence between unlike closed sections. Floater, *An O(h^{2n}) Hermite approximation for
conic sections*, 1997, section 2 and Theorem 3.1, gives geometric Hermite interpolants and a
Hausdorff bound for positive normalized middle weights up to 3; for degree three the subdivision
order is six. That bound controls the conversion of each section, not which points on two
sections should be joined. Repeated bisection of three initial parameter regions can explain a
knot pattern without determining that correspondence. No new pairing law is inferred from such
a pattern here. Source: <https://numerika.fmf.uni-lj.si/bib2html/e-documents/clanki/FloaterAnO(h2n)Hermiteapproximationforconicsections.pdf>.

Rational reparameterization does not supply a unique correspondence either. Hu, Wang and
Wang, *Piecewise Möbius Reparameterization of Rational Bézier Curves* (2018), section 2.1,
equation (3), changes weights by `w_i -> lambda^i w_i` under
`t = lambda u / (1 - u + lambda u)`, with `lambda > 0`, preserving degree and geometry.
Their section 2.2, equations (7), (9) and (12), uses curvature maxima as breakpoints,
normalized arc lengths as new breakpoints, and matching derivatives between pieces.
Those are choices for a single curve's speed, not a theorem selecting points to join
between unlike curves. In particular, a three-region knot pattern alone does not establish
that these regions are Möbius pieces. A numerical diagnostic fitting one such map through
three places of each region, then checking thirty other places, misses the observed conics
by 0.234 to 0.550 in model units across a five-set family (the maximum per section).
This is evidence against that three-piece hypothesis, not a bound on arbitrary piecewise
reparameterizations. The cross-section pairing remains unresolved; no new rule is installed.
Source: <https://www.jcad.cn/cn/article/pdf/preview/10.3724/SP.J.1089.2018.16764.pdf>.

The compatibility step must also be distinguished from choosing correspondence. Piegl and
Tiller, *The NURBS Book*, second edition (1997), section 10.3, pages 458-459, equation
(10.7), first assumes common section degree and knots, then interpolates their control
points, in homogeneous space for rational sections. Engleitner and Juttler, G+S Report 73
(2018), sections 2-4, permit incompatible section spaces through patchwork splines. Their
Theorem 3 proves existence in the chosen space from a blending construction under the
shadow-compatibility and smoothness-monotonicity assumptions; section 4.3 selects a solution
by constrained energy minimization. The sections' parameterizations remain given in both
constructions. These existence results therefore do not identify an unknown cross-section
pairing from an output knot pattern.
Sources: <https://home.zcu.cz/~bastl/GM1/the-nurbs-book.pdf>,
<https://www.gs.jku.at/pubs/NFNreport73.pdf>.


*Booleans.* Cell selection from one General Fuse with a ladder of fallbacks, each rung accepted only if
checked; a sheet tool divides rather than subtracts.

*Unproven.* Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean fallbacks, the sampling counts of the offset and
skin validity checks.


### Local focal witnesses for a sheet's normal sweep

`BRepOffsetExt_Focal` encloses non-periodic rational B-spline and Bezier supports directly from
stored poles, positive weights and knots, on each nonempty knot rectangle. The homogeneous first
and second derivative control nets use the derivative formula in Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing*, hyperbook edition (2009),
[section 1.4.2, equations (1.65)–(1.68)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
Their active-pole hulls enclose their values by nonnegativity and partition of unity
([section 1.4.1, (1.58)–(1.61)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node16.html)).
Each hull is intersected with interval Cox–de Boor evaluation on the current rectangle. Homogeneous
poles, location matrices, knot differences, quotient derivatives and fundamental forms all use outward-rounded
interval arithmetic; the calculation does not first convert or insert knots with ordinary floating
point arithmetic. This follows the arithmetic requirement in
[section 4.9.1](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node53.html).

Let I be the first fundamental form and II the second form against S_u cross S_v. Principal
curvatures are the generalized eigenvalues of (II,I); see the same authors,
[section 3.4, (3.41)–(3.43)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node30.html).
**Derived here:** the normal sweep S+s n, 0 <= s <= d (signed), has no focal point on a regular
parameter rectangle precisely when I-d II is positive definite. Indeed its eigenvalues in an
I-orthonormal frame are 1-d k_i. A negative diagonal entry or determinant proves that an eigenvalue
has passed zero somewhere along that normal segment. Sylvester's criterion proves the positive
case without subtracting nearly equal curvature roots. The signed curvature convention here is
opposite to that of the hyperbook's offset chapter; the focal condition is the same
([section 11.3.1, (11.18)–(11.25)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node224.html)).
The face orientation reverses d. A crossing witness must lie strictly inside the trimmed face;
it refuses thickening with kind `Curvature` and names that source face, before copying or building.

A regular verdict proves only local regularity within the support's knot rectangles. It does not
prove separation of distinct normal segments, nor continuity across knots. Periodic or other
supports, zero derivative denominators, singular parameter regions, exact equality and exhausted
work remain `Unknown`; a loose interval never causes a refusal. The existing construction remains
for unknowns, so complete general self-intersection exclusion is still open. The breadth-first
subdivision budgets (4096 rectangles, depth 24) are work choices, not geometric thresholds.

### A circular tube's parallel body after its hole closes

For a complete ring torus with core circle C of radius R and tube radius r<R, **derived here**:
T=C+B_r, so its outward distance-d parallel body is T+B_d=C+B_(r+d), using associativity of the
Minkowski sum and B_r+B_d=B_(r+d). Put a=r+d. At radial coordinate rho>=0 its meridian satisfies
(rho-R)^2+z^2 <= a^2, since the closest point on C has the same azimuth. When a>R, this is the
circle's disk clipped at the rotation axis. Revolving the circular arc through (R+a,0) between
(0,+/-sqrt((a-R)(a+R))) together with its axis segment constructs the boundary without the
self-crossing continuation. No surface fitting or new tolerance allowance is involved. The exact
hole-closure equality a=R remains on the existing route.

The construction applies to a valid positive-volume solid with one toroidal face. The source
face maps to the new skin; validity, positive volume and the existing independent clearance check
still apply. The independent volume oracle in the held-out family is the meridian's first moment:

V = 4 pi [ z0^3/3 + R (pi a^2/4 + R z0/2 + a^2 asin(R/a)/2) ],
where z0=sqrt((a-R)(a+R)).

This is a trimmed parallel **body**, whereas thickening retains the input sheet's normal sweep.
Trimming arbitrary spline offsets requires tracing all self-intersection components; the distinction
is described by Patrikalakis, Maekawa and Cho (2009),
[section 11.3.2, equations (11.28)–(11.29)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node225.html).

### Extending the approximation domain of a moved crease face

The kept-face offset approximates the exact normal offset on the UV rectangle of the
constructed face, including any extension needed to intersect adjacent moved supports.
It supplies that rectangle through `GeomAdaptor_Surface`, whose bounds select an evaluation
domain, rather than `Geom_RectangularTrimmedSurface`, which rejects bounds outside a bounded
basis surface. The approximation tolerance, degree and segment limits are unchanged.

Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and
Manufacturing* (2009), [section 11.3.1, equation (11.17)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node224.html)
defines the normal offset used here; [section 1.4.2, equation (1.62)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html)
expresses a spline as polynomial spans. **Derived here:** evaluating the end span beyond
its knot interval continues the same parametric support used by the extended edge curves;
it does not change the normal-offset definition or justify extrapolating a fitted result
outside its certified rectangle. The existing approximation-error, topology, clearance
and self-interference checks still apply to that rectangle. A generated family varies
scale, pose and signed distance on a valid near-planar C0 support, checking face lineage,
volume and independently evaluated normal-offset samples.

### Continuing both branches of a periodic seam

A seam has two curves on its face, although it has only one three-dimensional edge.
The two curves may differ by a period in either surface parameter. When the end faces
move, a finite seam pcurve can need the same end-span continuation as a tangent edge.
The kept-face construction now fits both branches independently, at its existing
one-tenth-tolerance interpolation bound, and gives the modifier the branch selected by
that edge use's orientation. Binding both curves by face identity alone would overwrite
one with the other. Branches that need no continuation keep their original route.
The solid validity, source-normal distance, self-interference and clearance checks remain.
A generated cylinder with finite spline seam pcurves checks the closed-form volume
pi (r+d)^2 (h+2d), at three scales and a rotated, translated pose, in both directions.

The continuation argument is **derived** from the existing parameter curves and the
normal-offset definition of Patrikalakis, Maekawa and Cho, *Shape Interrogation for
Computer Aided Design and Manufacturing* (2009), section 11.3.1, equation (11.17):
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node224.html>.
The representation contract is OCCT's `BRep_Builder::UpdateEdge(E,C1,C2,F,Tol)`
overload for a closed face, and `BRepTools_Modifier`'s separate `NewCurve2d` calls
for the two uses. This does not prove that arbitrary continuation is regular or
non-self-intersecting; those are still checked on the constructed body.

For an *exact* pipe P(s,t)=c(s)+r e(s,t), where e is the unit radial vector in the
normal plane of the regular spine, its normal is radial while the surface is regular.
Substituting into equation (11.17) gives c(s)+(r +/- d)e(s,t): its offset is a pipe
about the same spine. This substitution is **derived** from the pipe parameterization
and normal in the same authors' section 11.6.2, equations (11.111)-(11.115):
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node237.html>.
The local condition is that the shifted radius times the maximum spine curvature
is less than one; global separation also matters (section 11.6.3.4, Theorem 11.6.1):
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node242.html>.
A fitted tensor spline through circular sections is not thereby certified to be this
exact pipe. The kept-face operation has no retained spine/radius descriptor for a
general spline support. Inferring one from a nominal size and replacing the support
would change the input surface; this change does not do so.

### Collinear alignment: units and rational boundary rows

Collinear alignment uses exact degree elevation and the union of the two normalized
boundary knot vectors, with the larger multiplicity at each knot. It edits rows
directly and has no degree-five fitting restriction. Knot insertion preserves the
parameterized surface: Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer
Aided Design and Manufacturing*, hyperbook (2009), section 1.4.3, equations
(1.74)-(1.76), <https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node18.html>.
OCCT's `Geom_BSplineSurface::IncreaseDegree` preserves the geometry, while
`InsertUKnot` and `InsertVKnot` with `Add=false` raise a knot to the requested
multiplicity: <https://dev.opencascade.org/doc/refman/html/class_geom___b_spline_surface.html>.

**Derived here:** for unit row directions a and b with angle theta,
|a+b|=2 cos(theta/2). The symmetric fold check compares this dimensionless
quantity with 1e-9, a conditioning choice retained from the unit-scale algorithm;
physical short-row degeneracy is checked separately. Multiplying model lengths
cannot change this angular verdict.

Rational supports must have clamped boundary knots. Their rows carry weights as well as Euclidean poles. The moving boundary and
next row receive the reference's corresponding weights, after exact refinement;
deeper rows retain their weights. From the rational surface definition in the same
authors' section 1.5, equation (1.91),
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node20.html>, **derived here**:
write the boundary as S=A0/W0. Its cross derivative is
c (A1-S W1)/W0, where c is the positive endpoint degree/knot factor.
If the reference next-row weights are a constant multiple of its boundary weights,
then sum N_i w1_i P0_i = S W1. Consequently placing the moving next poles at
P0_i+r(P0_i-P1_i), with common r>0, gives the opposite proportional derivative.
This proves tangency for the constant-ratio option, including rational curved
boundaries. The proportional-weight condition is checked to 64 machine epsilons
relative (an arithmetic allowance); a nonproportional reference weight row with
that option remains unsupported. Without the option, pole lines align but neither
polynomial nor rational variable row ratios guarantee exact G1. Skew likewise
changes the target lines and is not covered by the zero-skew tangency argument.

At very large coordinates, independently evaluated exact isocurves and surface
points can differ by more than OCCT's default absolute confusion tolerance. The
Collinear face construction may carry the requested positional tolerance onto its
new face, edges and vertices if the default metadata does not validate. It then
rechecks BRep validity and the existing interior check. No poles, knots, weights,
curves or parameter curves are fitted or moved by this step, and no tolerance
larger than the caller's allowance is introduced. Already valid faces retain their
original metadata; the option is disabled for other consumers of the face builder.
The scale tests multiply both coordinates and the positional allowance by the
same factor, while checking boundary equality more tightly (1e-12 times scale)
and unit-normal agreement independently. Both planar and curved cases are tested
at scales 1e-3, 1 and 3e9, at rest and in a rotated, translated pose.

### Skin solid publication

Planar end caps are the only cap construction used by the smooth section skin.
A failed cap is a construction failure for a solid request; an open sheet request
continues to return its sheet. Before publication, closure is checked with
`BRepCheck_Shell::Closed`, orientation with `BRepLib::OrientClosedSolid`, and
validity with `BRepCheck_Analyzer`. Disconnected shells are not silently discarded.
These checks follow OPEN CASCADE, *OCCT Reference Manual* (8.0, 2026),
[`BRepCheck_Shell::Closed`](https://dev.opencascade.org/doc/refman/html/class_b_rep_check___shell.html)
and [`BRepLib::OrientClosedSolid`](https://dev.opencascade.org/doc/refman/html/class_b_rep_lib.html).
They establish topological closure and the kernel's validity conditions, not a
proof against geometric self-intersection. Nonplanar end-cap construction remains open.

### Rectangular sheet extension preserves trims

Extension requires a valid single-wire face whose boundary is its parameter
rectangle. `GeomAlgoExt::parameterRectangle` checks all boundary curves against
the four sides: line endpoints, or the complete positive-weight spline pole hull.
Area agreement cannot establish this condition, since a narrow notch can have
arbitrarily small area and finite depth. The pole-hull implication follows from
Patrikalakis, Maekawa and Cho, *Shape Interrogation for Computer Aided Design and
Manufacturing*, hyperbook (2009), [section 1.4.2, convex hull property](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
Applying that property to each constant coordinate is derived here. Unrecognized
boundary representations are not assumed rectangular.

### Coincidence of projected curves

Five common points and equal lengths do not imply equal curves. Plane projection
removes an edge only after establishing a common normalized spline basis (degree
elevation, knot insertion, and optionally reversal) with identical positive
weights. If every corresponding pole differs by at most the modelling tolerance,
the rational convex-combination formula bounds the displacement for every
parameter by that tolerance. The common parameter interval covers both edges,
so the bound is bidirectional. The pole-distance arithmetic is outward rounded.
This implication is derived here from Patrikalakis, Maekawa and Cho, *Shape
Interrogation for Computer Aided Design and Manufacturing*, hyperbook (2009),
[section 1.4.2](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
Representations whose correspondence cannot be established are both retained.
Degree elevation and knot insertion use OCCT's floating-point implementations;
the enclosure bounds the resulting common-basis representations, with the usual
rounding limitations of those representation conversions.

### Elementary focal distances

The offset Jacobian loses rank at a principal focal distance. Patrikalakis,
Maekawa and Cho, *Shape Interrogation for Computer Aided Design and Manufacturing*,
hyperbook (2009), [section 11.3.1, equations (11.24)-(11.25)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node224.html)
and [section 11.3.2, equation (11.28)](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node225.html),
give the offset curvature and singularity conditions. With this implementation's
normal convention, the factor is `1-d*k`. For a sphere or cylinder of radius R,
the inward focal distance is exactly R. These two elementary supports are classified
by comparing the signed physical distance to their stored radius, without the
relative allowance used for fitted blends.

For a cone with semi-angle alpha, `r(v)=R+v*sin(alpha)` and principal
curvatures are `0` and `-cos(alpha)/abs(r(v))`. Endpoint radii and a possible
zero crossing give the full focal range over the trimmed v interval, without
interior sampling. For a ring torus (`R>r`), the principal curvatures are
`-1/r` and `-cos(v)/(R+r*cos(v))`. Inward focal crossing is at r everywhere;
for an outward distance d the other crossing is `cos(v) <= -R/(r+d)`.
Endpoints and the contained multiples of pi give the cosine extrema over the
complete v interval. These reductions are derived here from the parametrizations
and the first and second fundamental forms. See Malakhaltsev and Vankerschaver,
*Sage Reference Manual: Differential Geometry of Parametrized Surfaces* (2010,
Sage 9.3 edition), [`principal_directions`, torus example](https://sporadic.stanford.edu/reference/riemannian_geometry/sage/geometry/riemannian_manifolds/parametrized_surface3d.html),
which gives these ring-torus curvatures with the same normal convention.
Arithmetic and trigonometric evaluation retain floating-point rounding limits.
Horn and spindle tori retain strict sampled classification.
A radius-10000 sphere moved inward by 9999.5 therefore remains a radius-0.5 sphere.
Classification of general fitted supports remains sampled with the existing
approximation allowance; it is not a proven curvature-error enclosure. A general
enclosure that separates uncertainty from a strict crossing remains open.

### Clearance coverage

Clearance carries separate coverage and resolution states. Empty input, failed
extrema, or exclusion of every region cannot become successful infinite clearance.
Distances are measured between whole trimmed regions and the original skin,
including boundaries; an interior lattice alone provides no lower bound.
Near openings, a region is either measured clear, or wholly inside the opening
exclusion, or subdivided. Work exhaustion is unresolved, never a measurement.

The exclusion bound is derived here: for a closed set A, the triangle inequality
gives `d(x,A) <= d(c,A)+|x-c|`. A containing box supplies c and a radius covering
the complete trimmed region. Thus `d(c,A)+radius < wall` excludes the entire
region. The radius is rounded upwards and the modelling confusion is added.
Trimming and distance computations still have OCCT's numerical tolerances; this
is a coverage check using those solvers, not interval certification of their
extrema. The existing distance allowance is unchanged.

A failed whole-region solve can instead use the derived lower bound
`d(x,skin) >= d(c,skin)-radius`, subdividing until that bound resolves the
requested clearance. This uses the same triangle inequality. A failed point
solve or exhausted subdivision stays unresolved.

Distance queries use independent copies of each trimmed surface. Spatial edge
curves inherited from an import can deviate from that surface within their
stored tolerance. Measuring those curves would measure their representation
error as wall thickness. Rebuild query curves from each face's own pcurves and
support, retaining the published surfaces and trims. `GeomLib::BuildCurve3d`
returns a maximum approximation deviation; unresolved or excessive deviation
fails the query. For query curve error at most e on either shape, the triangle
inequality gives `|d(A,B)-d(A',B')| <= 2e` (derived here). Subtract this amount
from the numerical distance; add it to an opening exclusion upper bound.
Here e is one tenth of the existing modelling confusion. This remains a
numerical approximation/error estimate, not a Bernstein-certified enclosure
of the approximation. OCCT's *Boolean Operations* specification, section
"Interferences", describes the stored-tolerance model; a stored edge tolerance
is not a justified allowance to reduce the requested wall thickness.

### Walls and opening rims

A shell's clearance domain is the inward descendants of the retained input
faces. Descendants reached only from removed input faces close the opening;
they are rims, not offset walls. This follows the construction in Open Cascade,
*Open CASCADE Technology Reference Manual* (2026, version 8.0.1),
[`BRepOffsetAPI_MakeThickSolid`, Detailed Description and `MakeThickSolidByJoin`](https://occt3d.com/dev/doc/refman/html/class_b_rep_offset_a_p_i___make_thick_solid.html):
the first shell contains the input faces except the removed faces; the second
is its offset. Track this distinction by input-face identity through construction,
not by surface type, coordinates, a distance threshold or a hash. A descendant
shared with any retained input face is still measured. Opening descendants
remain in the reported history. An empty retained domain still fails coverage.

For example, an opening in a rectangular pocket can produce an annular rim
slightly beyond the removed face, with an outer trim at the requested offset
and an inner trim touching the original pocket walls. Independent whole-face
distances give zero against those walls. No positive lower distance bound
exists on this rim. A uniform interior lattice can miss that zero.
For a plane separation h and an in-plane distance r from the rectangular
opening, `d_open^2 = r^2 + h^2`, while the adjacent trimmed wall can be at r.
Thus a spherical opening exclusion of radius t retains distances down to
`sqrt(t^2-h^2)`, not t (derived here by Pythagoras, where h < t).
Subdivision cannot repair this difference of domains.

The regression test uses a pocketed body with a bore, a round and a chamfer,
at three scales, two poses and three wall thicknesses. Separate trimmed-face
extrema measure all 17 retained wall descendants at the requested thickness;
the opening descendant has minimum zero and remains in the history. A separate
thin retained face fails. The numerical solver and its tolerance are documented
by Open Cascade, *Open CASCADE Technology Reference Manual* (2026, version
8.0.1), [`BRepExtrema_DistShapeShape`, Detailed Description and `Perform`](https://occt3d.com/dev/doc/refman/html/class_b_rep_extrema___dist_shape_shape.html).
These are numerical distance checks, not interval proofs of extrema.

### Regularity of network surfaces

A lattice sized from the number of knots does not cover short knot intervals.
Rejection witnesses now include every actual knot interval. Acceptance additionally
uses the outward-rounded Cox-de Boor derivative enclosure on every knot rectangle
implemented by `focalSweep` at distance zero. The derived identity
`det(I) = |Su cross Sv|^2` permits the cross-product enclosure to establish
regularity directly without dependency loss in `E*G-F*F`. A zero-free interval
normal proves local immersion; unresolved cells do not receive that verdict.
The basis and convex-hull bounds are those of Patrikalakis, Maekawa and Cho,
*Shape Interrogation for Computer Aided Design and Manufacturing*, hyperbook
(2009), [sections 1.4.2 and 1.4.3](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html).
Periodic finite covers use the stored extended flat knot sequence and cyclic pole
copies, matching OCCT's BSplSLib evaluation without fitting. This proves local
immersion only. A separate global injectivity certificate remains open; the
argument analyzer can flag valid periodic boundary seams and does not establish
self-intersection freedom within a single face. It is not an acceptance proof here.

To reduce interval dependency for rational surfaces, translate homogeneous
coordinates by a constant local anchor: `(X,W) -> (X-aW,W)`. This changes the
surface to `S-a`; every first and second derivative is unchanged (derived here).
All subtraction and multiplication are outward rounded. At distance zero only
first derivatives are needed, so second derivative nets are not evaluated.
Local-immersion queries have no work-count cutoff; the existing finite
subdivision precision limit still reports unresolved, never regular.

### Neighbour-bounded planar push-pull

Zou and Feng, *Push-Pull Direct Modeling of Solid CAD Models*, preprint (2018),
[section 5.2, Algorithm 3](https://arxiv.org/pdf/1811.04036), construct auxiliary
volumes bounded by the neighbours' extended underlying surfaces and the moved
face's start and end supports. A prism is unsuitable for a slanted neighbour.

The convex planar specialization here is derived: write the body as an
intersection of outward support inequalities `n_i dot x <= c_i`. Translating the
selected faces by v changes only their constants to `c_i+n_i dot v`. The result
is their new halfspace intersection. This is exactly the neighbour-bounded volume
update; active constraints may appear or disappear without a special topology
rule. Moving constants leaves the recession cone unchanged, hence a bounded
original intersection stays bounded whenever the target is nonempty. Triples of
independent support planes contain every vertex, providing the temporary box.
No artificial box face may survive in the result. Each surviving face is tracked
through Boolean history back to its defining original support.

Before using this construction, both Boolean differences between the original
body and its reconstructed support intersection must contain no solid. This
rejects nonconvex input. Final topology must pass the normal validity check.
The construction is attempted only after the existing sweep rejects a slanted
neighbour; previously successful sweep geometry is unchanged. Curved supports and
nonconvex neighbour-bounded auxiliary volumes remain open.

*Adjacent refinement (derived).* A rejected adjacent fit marks angular and
complete-form curvature misses as well as positional misses. Between-row side
checks call the same measurement and easing alternatives as acceptance. A fit
already accepted retains its previous refinement decisions. This removes a
position-only stop condition; finite refinement still does not prove existence
or convergence for arbitrary boundary jets.

Endpoint jet information is retained as a hint after a fit has already failed its
result checks. It explains the failed fit and never rejects a delivered fit or
acts as a feasibility condition for positional-tolerance guides.

### Owned point witnesses near curved trims

A fixed-resolution trim polygon can misclassify points near a curved rim. An
interior-distance witness now uses parity of boundary arcs certified to avoid
the query point. Positive-weight Bezier hulls enclose each arc in a convex tube
about its chord. If that tube excludes the point, straightening the arc to its
chord cannot cross the point and preserves crossing parity. Subdivide only
uncertain arcs; retain uncertainty once hull deviation reaches the existing
numerical enclosure guard. For a rectangle witness, additionally require that
no boundary enclosure intersects the rectangle. Surface construction, result
bounds and the ordinary surface-interior checks remain unchanged.

The hull property is Patrikalakis, Maekawa and Cho, *Shape Interrogation for
Computer Aided Design and Manufacturing* (2009), section 1.4.2. Polygon crossing
parity and its relation to winding number are given by Hormann and Agathos,
*The Point in Polygon Problem for Arbitrary Polygons* (2001), sections 1-2,
equations (1)-(3), [author manuscript](https://www.inf.usi.ch/hormann/papers/Hormann.2001.TPI.pdf).
The adaptive arc-to-chord restriction is derived here. This is not interval
arithmetic: it retains the established floating-point enclosure guard and treats
unresolved curves or near-boundary points conservatively.

Adjacent owned pcurves can have different representations of their common
vertex. A parity chain includes connecting segments only when each gap fits
within the unchanged enclosure guard. For larger gaps, point ownership comes
from `BRepClass_FaceClassifier` on the delivered face, rather than an invented
closing curve. Before propagating that ownership across a nonzero rectangle,
exclude every boundary-arc tube and the entire coordinate box of every endpoint
gap, enlarged only by the existing guard. Unsupported pcurves cannot supply a
cell certificate. A remote endpoint gap therefore does not invalidate an
otherwise resolved local cell. Points classified on the rim remain unresolved.

This local propagation is derived: an open connected rectangle disjoint from
the delivered boundary lies in one component of its complement, so its ownership
is constant and one authoritative interior point suffices. It assumes the face's
own ownership semantics; it does not certify a malformed wire's missing geometry
or assert that any arbitrary-gap closure is the delivered trim. The test family
varies endpoint gaps, UV scale and rigid pose, and independently checks the seed
with the face classifier. The ordinary interior-curvature classifier, all
construction arithmetic and all result tolerances are unchanged.
