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
criterion above. v2 now makes all four at 70 and the two at 75 as the balls' envelope (every point r
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
linear, falling, peak and valley, Round and Conic rho 0.3) v2's lie within 1e-12 to 1.1e-6 mm both ways;
round a closed rim the alternating law within 8e-6, but a law of three uneven places (1, 2, 1.5) 0.04 off,
the reference modeller's closed interpolation being neither the periodic nor the natural spline (open). Native
(`BRepBlendExt_Test`): rounds on a box edge (linear, steep, peak, valley) lie on their circles to 1e-6
and remove (1 - pi / 4) times the integral of r^2 to 1e-9 of the body, at rest, turned and at scale 1e-3
and 1e3; a cylinder's rim (periodic) and, held out, the ring of a sphere on a plane on their circles to
1e-6; chamfers of a growing setback, also at 30 degrees off a reference face and within a range, remove
the integral of their sections to 1e-9. G2 in the application (the audit's variable
and range rows): seams G1 to 0.01 degrees and G2 to 1e-4 / r.

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
than the ball; v2 makes it (the G2 section cut by the bottom face, its seams G2 to 1e-14 / r).
A hole tangent to a face beside the edge pinches the body: the hole meets the face in a line along
which four faces meet. No blend passes it and v2 refuses it by name, as the reference modeller does; with a gap of
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
between nets, up to the finest).
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

*Size.* Lengths are shares of the loop's size or of the tolerance, with one exception by design: where a
fill goes over into planes and other faces at once (or into a cone alone) its rim and its edges' parameter
curves are made to a tenth of `Precision::Confusion()`, OpenCascade's absolute least distance, so that its
edges keep that as their tolerance and STEP gives the joint back (a G1 cap over a circle on a cone, tolerance
1e-4 in step with the size, is then made on 11 x 11 poles at 1e-3 and on 17 x 17 at 1 and 1e3). A fitted net moved onto its curves ends
its ladder when the bounded rim correction leaves its own places further off than the tolerance (the
correction's own residual, no constant). The fills of the tests, with every length and the tolerance
scaled by 1e3, are made or refused alike with the same nets; at 1e-3 so are a disc, a cap, an open loop of
arcs, the dish at G1 and the elbow. Open: at 1e-3 the dish at G0 (its tolerance then
`Precision::Confusion()`) is refused as folding over, the corner of three rounds (its tolerance below it)
is made, the square on rounds misses G2 by a share of 0.109 by its quad and is trimmed instead, the
pentagon and the split side take other nets; the teddy takes another net at 1e-3 and is refused at 1e3. A
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
tolerance). *Sharp edges:* the new edge is the crossing of the two offset surfaces, solved in the plane
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

*Booleans.* Cell selection from one General Fuse with a ladder of fallbacks, each rung accepted only if
checked; a sheet tool divides rather than subtracts.

*Unproven.* Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean fallbacks, the sampling counts of the offset and
skin validity checks.

