<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# Limits and open work

What the extensions do not do yet, with the cause of each. The top-level [README](../../README.md) lists the short
form; the derivations and the unproven choices are in the [operations reference](operations.md).

## Blend solver

An edge whose faces meet at a steady angle is blended at any opening from 5 to 175 degrees, convex or concave.
Convexity is proven along the whole edge in Bernstein form for spline and Bezier surfaces, extrusions of spline
curves, planes, cylinders, cones, spheres and tori; revolutions and offsets are read at the places k/16 of the
edge, which is sampling. An edge with convex and concave parts is refused by name.

Where the angle turns along the edge so far that the balls' centre path bends more tightly than the ball (the
rim of a hole leaning at a with `r >= (R + r) cos a`: from about 70 degrees), the rolling-ball surface crosses
itself. A contact that crosses itself is now found in its face's parameters on any surface and cut at the
crossing: a round's loop is patched by its ball, and a chamfer's two sheets meet in their crease (measured on
the slanted rims of a plate, leaning 70 to 85 degrees: all 16 cases chamfered).

Where the spine bends more tightly than the ball, a round's sheet crosses itself between its contacts and is
trimmed along that crossing to the balls' union (see [the section](operations/BRepBlendExt.md#where-a-round-crosses-itself-between-its-contacts));
other laws are still refused there by name. A fixed-contact conic replacement on a cylindrical
cavity and transverse plane retains the envelope outside its crossing sector; continuous spline and embedding
certificates for it remain open (see [the conic-window section](operations/BRepBlendExt.md#fixed-contact-conic-windows)).

Isolated shared-edge events between joined contact strips are partitioned; overlapping or simultaneous boundary
events need a more general contact arrangement (see [contact
ownership](operations/BRepBlendExt.md#contact-ownership-on-a-joined-blend-sheet)). A blend ends on a face
reached past the vertex only where that face's surface also meets the other contact.

Some configurations where several blends meet are refused although a valid result exists: trimming a blend group
needs an owner for every retained region and contact before chain ends are cut, and the centre path is
parametrised by arc length along the traced offset spine only for an open constant-distance rail whose own
stations fold.

The later construction passes (extended fits, refined contacts, offset spines) each run only where the pass
before read the input they change; where one refuses and a later pass measured the obstruction on its sheet, the
refusal names that obstruction.

Corners where four or more blends meet: conic sections are built, split into exact rational patches of degree 6
that are G1 to the sheets, where every face at the corner is a plane (all 144 members of an uneven family build;
the corner bends more than ten times as tightly as its blends on 22 of them, all on one steep pentagon).

Curvature-continuous and chamfer corners of four or more blends, and conic ones with a face that is not a plane,
are refused (see [split corners](operations/BRepBlendExt.md#split-corners-of-more-than-three-conic-blends)).
Across a knot where a spline support is only C2 the curvature-continuous section is smoothed in a window derived
from leading-order secants, not certified suprema; where the spans do not meet at the knot the exact law is
kept.

A curvature-continuous blend must meet the current seam and section-law bounds. A limited blend is built
directly from the blend's faces only under the premise stated in [the variable-size
section](operations/BRepBlendExt.md#a-size-that-changes-along-the-edge) (one edge, both limits inside it, every
other changed edge outside the slab, checked by sampling); otherwise the Boolean isolation stays.

## Section law for curvature continuity

`BRepBlendExt_Policy` supplies the conic section only. A published bound for a curvature-continuous conic or
quintic section as a function of the support's curvature is missing.

## Fill

The interior acceptance gate and the flat-support curvature limit have no derivation. Across two rims of an open
boundary the gate's reference is the Hermite bridge between them, which is reasoned; the factor it allows over
it is a judgement.

The ladder of nets stops where two nets in a row agree; its rungs are not nested spline spaces, so on some loops
the distance between nets does not shrink (a plateau, or nets that never settle) and the stop is then a
judgement, and a guide's miss can rise from one rung to the next.

At a corner whose faces bend differently the bend asked eases to the one a surface can have over one span of the
coarsest net, a judgement: G2 is met on the middle of each side only, the corner's miss is reported, and the
coarse nets of such loops can be lopsided.

Exact G2 at a regular corner whose supports have different second fundamental forms is impossible, since
continuity would require the fill's one form there to equal both (the [corner compatibility
proof](operations/BRepPatchExt.md#g2-at-a-corner)); this includes a bending strip meeting a flat strip and
rounds of unequal forms, and neither local refinement nor several internally G2 faces removes the obstruction.

The eased fill, with its reported corner miss, relaxes that request; some non-four-sided loops are still refused
because the finite net misses the flat-support limit beyond the ease, which the proof does not show impossible.
Fills through three profiles with curvature continuity can miss the curvature-share limit. Open boundaries are
chained by `BRepPatchExt_Chain` (nearest-end growth, gaps closed with straight segments); which chain a user
means is not decided in general.

A loop that is a height over no plane is filled over a faired disc laid flat conformally (Boundary First
Flattening); the net over that chart is a tensor net, so a loop that needs fine spans in one place needs them
right across. Its 640-pole forecast is a computation policy based on osculating-circle derivatives, not an
approximation bound or a necessary pole count. Open: complete curvature enforcement, converged geometric
feasibility and support-derived collar initialization.

Not supported: curvature-continuous fills between plates meeting at three or more edges at once, strips with
bending flows, flow-direction modes on steep supports, and profile-curve panels.

## Offsets and thickening

Strict local focal crossings on non-periodic rational spline sheets are refused by interval witnesses; periodic
supports, inconclusive bounds and global collisions are unproven. A complete circular tube's outward parallel
body closes its hole by revolving its clipped meridian. Offsetting solids whose blends meet the wall thickness
exactly has no exactness proof. In the face offset the largest angle between the faces along a tangent edge is
found by sampling, not bounded.

A separated tangent junction can add a ruled contact strip after the fixed-incidence route refuses; general
continuation of approximately tangent contacts is open, and a strip needs separate sharp-edge traces and a valid
closed arrangement. Endpoint gaps are checked independently of result tolerances, both seam pcurves included.

See [the endpoint certificate limits](operations/BRepOffsetExt.md#face-offset-endpoint-certificate-limits) and
[separated tangent junctions](operations/BRepOffsetExt.md#separated-tangent-junctions).

## Speed

The ratios against stock followed by the same check are in [Results](../../README.md#results) and
[benchmarks.md](benchmarks.md), each with its cause. The test that a fitted blend does not cross itself is
evaluated on every fitted blend; a curvature bound that proves the sections cannot run back would let it be
skipped (not done).

## Tests

Operations without a test of their own yet: deformation of whole bodies (`deformBody`, `deformSupports`),
Boolean imprints along faces' normals, and the curvature-continuous fills between plates meeting at three or
more edges (not yet supported). Of the tests that exist, the ones for regions, thickening, bridging and
aligning were each shown to fail with the code or the continuity order broken; the others have not been
mutation-checked. The test suite is described under [Tests](../../README.md#tests).

## Conventions

Formatting is OCCT's (`clang-format` 18.1.8 with the repository's file), headers carry the repository's licence
text, comments use `//!`, and files are named `Package_Class`. Not done: the packages use free functions in
namespace `BRepExt`, the handle interface and `std::runtime_error` failures where OCCT style would be classes
with `IsDone` / `Shape`, `Standard_Failure` and `Message_ProgressRange`. Public free functions carry
`Standard_EXPORT`; a Windows DLL build has not been tried.

## Policy

The choices listed under [Section policy](../../README.md#section-policy) are not mathematics and are not
included.
