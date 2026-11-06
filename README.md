<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# OCCT extensions

**A modelling layer for [Open CASCADE Technology](https://dev.opencascade.org) 8: blends that are exact where stock's are approximate, fills and skins with real continuity control, offsets that survive vanishing faces, and twenty-three fixes to the algorithms underneath.** Every result is checked against a definition or a closed form, and every operation checks its own answer before returning it.

The repository is OCCT with its history, followed by commits that add

* **`TKAlgoExt`**, about 60 000 lines in nine packages: rolling-ball fillets and chamfers, boundary-first Patch fills, Gordon and guided skins, framed sweeps, offsets and shells, G0 to G3 bridge curves, surface alignment and squaring, curve and surface deformation, and cell-based Booleans;
* **numerical utilities** in `TKGeomAlgo` and `TKMesh`: interval enclosures, Bernstein products, constrained tensor-spline fitting, certified knot compaction, trim integration, surface jets;
* **twenty-three corrections** to existing OCCT algorithms, one commit each, with the derivation and limits of each. One of them makes face classification 30 times faster on the benchmark plate.

The base is OCCT `b8f597c6` (30 July 2026). Nothing here has been proposed upstream yet; it is organised so that it can be proposed in pieces ([Upstreaming](#upstreaming)).

## What you get that stock OCCT does not give you

| Capability | Stock OCCT | This tree | Evidence |
| --- | --- | --- | --- |
| **Rolling-ball fillets** on curved supports, chains of rounds and corners | `BRepFilletAPI_MakeFillet` leaves a defect of up to 7e-3 where a round runs into earlier rounds and 3e-3 on crossed cylinders | Blends exact to the sampling (defect 1e-8 to 1e-14); radius laws along an edge, ranges, stations, chamfers, corners of three or more blends, junctions of two curved surfaces | [Results](#results), `BRepBlendExt`, `BRepFilletExt` |
| **Patch**: fill a loop of edges with G0, G1 or G2 to the neighbouring faces | `BRepFill_Filling` / `GeomPlate` plate fit: 0.11 to 0.45 off the defining surface on four cubic sides | Fill 1e-12 from it; n-sided loops, open boundaries chained and closed, interior guides and points; loops that are a height over no plane (rims whose faces sweep a half circle or more) filled over a conformal chart of their own; the same net at a thousandth and a thousand times the size; a G1 request that is geometrically impossible is refused with the reason | [operations.md](dox/extensions/operations.md#brepatchext-boundary-first-fill) |
| **Smooth skins** through sections, and Gordon skins through networks of profiles and guides | `BRepOffsetAPI_ThruSections` bulges to half-width 19.3 through six circles of radius at most 10 | Stays within the sections (10.0); closed rings joined seam to seam; guided skins that contain their guides | [benchmarks.md](dox/extensions/benchmarks.md) |
| **Framed sweeps and drafts** | `MakePipeShell` twists a profile unless the path starts square to it; `DraftAngle` takes planes, cylinders and cones | Rotation-minimising frames (double reflection, fourth-order error); draft of any wall by its isocline surface | [operations.md](dox/extensions/operations.md) |
| **Hollow, shell, offset, thicken, push-pull, defeature** | An offset that reaches a fillet radius returns an invalid or wrong solid (volume off by 3.3e4 on the benchmark) | Faces the distance would shrink to nothing are removed first; exact to 1e-11 on the benchmark; a face offset that keeps every face, sharp edges where the moved faces cross, tangent edges held to the input's own tolerance, each result certified by distances and loadable from a saved file; sheets with C0 creases are divided before they are thickened; several features removed in one call | [benchmarks.md](dox/extensions/benchmarks.md) |
| **Bridge curves of G0 to G3** between curves or from a point and a direction | no direct counterpart | Continuity from the end conditions of the geometric-continuity constraints (Hohmeyer and Barsky 1989), tested to 1e-12 on a circle and a helix; cusps refused | [references.md](dox/extensions/references.md) |
| **Surface bridging, alignment, squaring, deformation** | no direct counterpart | G1 bridges and alignments whose normals agree with the supports along the shared edge; four-sided surfaces over boundaries with supports; carrying curves and bodies from one face to another in UVN | `BRepSurfaceExt`, `BRepCurveExt` |
| **Curves**: offset at any scale, outline, regions, corner rounding of wires, tangent circles, spirals | `BRepOffsetAPI_MakeOffset` depends on the curve's size (absolute 1e-7 lengths); `HLRBRep` has documented limits on tangent and edge-on faces | Offsets that behave at any scale, exact silhouette events, all bounded cells of overlapping and nested loops | `BRepCurveExt` |
| **Booleans that report what they did** | a single run | A ladder of fallbacks, each accepted only if its result is checked; imprints of curves along a direction or normals | `BRepBooleanExt` |
| **Certified numerics** | no direct counterpart | Outward-rounded interval enclosures, a displacement bound for knot removal over the whole rectangle, adaptive trim integration | [utilities.md](dox/extensions/utilities.md) |

Where this tree is slower than stock plus the same check, or does not do something, [Results](#results) and [Open work](#open-work) say so with the cause.

## Results

![CPU time per call: stock OCCT, stock followed by a validity check, and TKAlgoExt](dox/extensions/operation-benchmark.svg)

Median CPU time per call of ten runs, the median of three such runs. Every call of `TKAlgoExt` checks its own result before returning it, so the
fair reference is the middle bar, stock followed by the same check. Hatched bars are calls with no valid or a
wrong result. Selected rows (milliseconds; full table, method and caveats in
[benchmarks.md](dox/extensions/benchmarks.md)):

| Operation | Stock | Stock + check | TKAlgoExt | Quality |
| --- | ---: | ---: | ---: | --- |
| Box, all 12 edges rounded | 2.6 | 3.5 | 1.8 | both exact to 3e-11 (Steiner) |
| Sphere R14 on a box, junction rounded | 1.02 | 1.23 | 0.48 | ball defect: stock 4e-6, TKAlgoExt 1e-14 |
| Round blended into earlier rounds | 5.8 | 7.3 | 8.6 | ball defect: stock 7e-3, TKAlgoExt 2e-9 |
| Crossed cylinders R10 / R6, r = 1 | 10.1 | 12.2 | 11.6 | ball defect: stock 2e-4, TKAlgoExt 2e-8 |
| Fill of four cubic sides | 2.7 | | 2.0 | stock 0.11 off the defining surface, TKAlgoExt 1e-12 |
| Hollow of a filleted box | 2.6 | | 6.7 | stock volume wrong by 3.3e4; TKAlgoExt exact to 1e-11 |
| Smooth skin through six circles | 1.07 | 1.26 | 0.99 | stock bulges to half-width 19.3; sections have radius at most 10 |

It is not uniformly faster: a single circular fillet and the crossed cylinders are level (1.05 and 0.95 times stock plus
check), the round blended into earlier rounds is 1.19 times, the curve offset 1.5 times, the smooth skin through square,
circle and square 4.3 times, and the hollow 2.6 to 3.2 times slower while returning a correct body. The rounds pay for a
test that a blend does not cross itself, which steep rims need. [benchmarks.md](dox/extensions/benchmarks.md) gives the
cause of each.

## Gallery

Results of the operations, rendered from the meshes of example models with three.js.

| | |
| --- | --- |
| ![Box with every edge rounded](dox/extensions/images/fillet-corners.png) | ![Rounded body hollowed to a wall, cut open](dox/extensions/images/hollow-cutaway.png) |
| Twelve edge blends meeting at eight corners (`BRepBlendExt`, `BRepFilletExt`). | The same body hollowed to a wall of 2 and cut in half (`BRepOffsetExt`). |
| ![A five-sided fill between free-form faces](dox/extensions/images/mixed-freeform-5.png) | ![A fill at a corner of unequal, set-back rounds](dox/extensions/images/mixed-setback-unequal-rounds.png) |
| A five-sided fill between free-form faces (`BRepPatchExt`; the neighbouring faces are the grey strips). | A fill at a corner of unequal, set-back rounds; the orange face joins the grey rounds smoothly. |

| | |
| --- | --- |
| ![A small round blended into large rounds](dox/extensions/images/fillet-reblend.png) | ![Three rounds meeting at a box corner](dox/extensions/images/fillet-three-edge-corner.png) |
| Four small rounds on the uprights run into the four large rounds of the top face (rounds blended into existing rounds). | Three rounds of one radius meeting at a single corner. |
| ![A fill at a corner of three unequal rounds](dox/extensions/images/mixed-unequal-rounds.png) | ![A fill between a plane, a cylinder and a sphere](dox/extensions/images/mixed-plane-cylinder-sphere.png) |
| A fill at the corner of three rounds of different radii; the orange face is tangent to all three. | A fill whose boundary lies on a plane, a cylinder and a sphere. |
| ![A three-sided free-form fill](dox/extensions/images/mixed-freeform-3.png) | ![A six-sided free-form fill](dox/extensions/images/mixed-freeform-6.png) |
| A three-sided fill between free-form faces. | A six-sided fill between free-form faces. |
| ![A four-sided fill between spline faces](dox/extensions/images/mixed-bspline-4.png) | ![A fill between a round and a flat face](dox/extensions/images/mixed-fillet-flat.png) |
| A four-sided fill between spline faces. | A fill between a round and a flat face. |

| | |
| --- | --- |
| ![Two crossing cylinders with the intersection rounded](dox/extensions/images/junction-crossed-cylinders.png) | ![A sphere on a box with the junction rounded](dox/extensions/images/junction-sphere-on-box.png) |
| Two crossing cylinders (R10, R6) with their intersection curves rounded by 2.5: the blend follows a space curve where two curved surfaces meet. | A sphere (R14) standing on a box with the circular junction rounded by 6. Stock and this branch agree on the volume to 1e-8 of it. |

| | |
| --- | --- |
| ![A fill between three tube ends at a Y junction](dox/extensions/images/model-y-pipe.png) | ![A fill bounded by the edges of six plates](dox/extensions/images/model-mouse-network.png) |
| A fill between three tube ends at a Y junction, curvature-continuous with the three grey sheets. | A fill bounded by the edges of six plates round a three-way junction. |
| ![An open fill between three strips](dox/extensions/images/model-form-finding.png) | ![A cap between cylindrical surfaces](dox/extensions/images/model-cylinder-cap.png) |
| An open fill between three strips (the boundary is not a closed loop), shaped by guide curves. | A cap between two cylindrical sheets, the whole model turned in space. |
| ![A hull bottom between sheer lines](dox/extensions/images/model-hull-bottom.png) | ![A domed top on its boundary curves](dox/extensions/images/model-mouse-top.png) |
| A hull bottom between sheer lines, turned in space. | A domed top on its boundary curves, turned in space. |
| ![A panel joined to an existing neighbour](dox/extensions/images/model-bumper-panel.png) | |
| A panel built on an existing neighbouring panel. | |

The last block shows parameterised models of larger shapes with several faces around one fill, open boundaries
and guide curves, each tested at two poses against a closed-form or symmetry oracle. The images illustrate what
is supported; the numbers are in the benchmark and the tests.

## What is added

| Part | Where | What it is |
| --- | --- | --- |
| Corrections | one commit each, [corrections.md](dox/extensions/corrections.md) | Fixes to existing algorithms: quadratic or unbounded loops, wrong ranges, uninitialised reads, missing cancellation. 17 build against current `master` (one needs another of them), four later ones are listed with them. |
| Numerical utilities | `GeomAlgoExt` in `TKGeomAlgo` | Interval enclosures, Bernstein products, spline collocation, constrained tensor-spline fit, rim fit, certified knot compaction, surface jets |
| Patch chart and surveys | `BRepAlgoExt` in `TKMesh` | Cotangent-weight chart of a boundary loop, interior fairness survey, trim integration, read-back |
| Operations | `TKAlgoExt`, nine packages, about 60 000 lines | Listed below |

| Package | Operations | Lines |
| --- | --- | ---: |
| `BRepBlendExt` | Rolling-ball blend solver: ball on two supports, spine marching, section fit, attachment to the body | 15 800 |
| `BRepFilletExt` | Constant, variable, range and station fillets, chamfer, contour and prism support, blend validity checks | 4 900 |
| `BRepPatchExt` | Patch: one face inside a loop of edges (chained, with gaps closed, if the boundary does not close), with G0/G1/G2 to neighbouring faces and interior guides, over the mean plane or over a chart of its own | 8 100 |
| `BRepSurfaceExt` | Surface alignment, bridging, constraint jets, deformation, face support and rebuild, face separation proofs | 7 900 |
| `BRepCurveExt` | Planar and 3D curve offset, outline, tangent circles, corner rounding of wires, G0 to G3 bridge curves, projection, rebuild, refinement, regions, spiral | 6 600 |
| `BRepOffsetExt` | Solid offset, hollow, shell, thicken, push-pull, face offset that keeps every face, defeaturing, sewing history | 5 000 |
| `BRepSkinExt` | Gordon and guided skins, framed sweeps, draft | 4 500 |
| `BRepBooleanExt` | Cell-based Booleans with a ladder of fallbacks, checked solids, imprints | 3 900 |
| `BRepToolsExt` | Properties, bounding boxes, dimensions, fairness, hidden-line, failure hints, shared checks and history | 2 800 |

Each operation is described, with its construction, mathematics, sources, unproven choices and stock
counterpart, in [operations.md](dox/extensions/operations.md). The utilities and their error budgets are in
[utilities.md](dox/extensions/utilities.md), the literature in [references.md](dox/extensions/references.md).
"Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the
source next to their use.

## Using the operations

The operations keep shapes under 32-bit handles and report failures by exception. The interface between an
embedding program and the operations is four declarations in `BRepToolsExt_Session.hxx`, which the embedding
program defines:

* `store(OcctKernel&, const TopoDS_Shape&)` registers a shape and returns its handle;
* `get(const OcctKernel&, handle)` returns the shape behind a handle;
* `release(OcctKernel&, handle)` drops a temporary;
* `sessionCancelled()` is polled for cancellation.

`OcctKernel` is declared and never defined here; it is only passed back. A program that keeps shapes in a
`std::vector` needs twenty lines (`adm/scripts/extensions/operation-benchmark.cxx` and the test support
`TKAlgoExt/GTests/BRepExtTest_Support.cxx` are two). Failures are thrown as `BRepExt::Failure`, which derives
from `std::runtime_error` and carries a stable reason code.

### Section policy

What follows from the mathematics is in this tree. What an application chooses is not: the blend section law
for curvature continuity, how a bridge curve chooses tensions it was not given, which control vertices an
alignment keeps, the shape of a sheet extension, and the mapping of a curve deformation. For curvature
continuity, `BRepBlendExt_Policy.hxx` has a conic-section default and refuses curvature-continuous requests
until a caller supplies a law through `BREPBLENDEXT_POLICY_HEADER` and `BREPBLENDEXT_POLICY_CLASS`
([Building](#building)). The derived parts that those choices sit on are here: the G0 to G3 bridge
construction `BRepCurveExt_Bridge.hxx`, which takes the tensions as arguments, and planar corner rounding
`BRepExt::roundWireCorners`.

## Building

The tree is an ordinary OCCT CMake project. `CMakePresets.json` holds two configurations of the modules the
extensions need (no Visualization or Draw, static Release libraries, no FreeType or RapidJSON):

```sh
cmake --preset native && cmake --build --preset native            # build machine
source /path/to/emsdk/emsdk_env.sh                                 # Emscripten 5.0.3
cmake --preset wasm && cmake --build --preset wasm                 # WebAssembly
```

The `wasm` preset sets the toolchain file from `$EMSDK` and compiles with
`-fwasm-exceptions -O3 -msimd128 -DIGNORE_NO_ATOMICS=1 -DOCCT_NO_PLUGINS`. Applications may also
`add_subdirectory()` the tree with the same variables set, and then define the session functions above in
their own objects.

Two options belong to the extension packages (`adm/cmake/occt_extensions.cmake`):

* `BUILD_EXT_KEEP_EXCEPTIONS` (on by default) compiles `TKAlgoExt`, `GeomAlgoExt` and `BRepAlgoExt` with C++
  exceptions and assertions on whatever OCCT's own flags say (`-UNo_Exception -UNDEBUG`), because the operations
  report failures as exceptions and rely on their checks.
* `BREPBLENDEXT_POLICY_HEADER` with `BREPBLENDEXT_POLICY_CLASS` name a header and a class that replace the
  built-in section policy of the blend solver (force-included into the `BRepBlendExt` sources).

## Tests

`-DBUILD_GTEST=ON` adds the target `OpenCascadeGTest`, which runs the tests of every built toolkit; 116 tests
belong to this tree (92 for the operations, 16 for the utilities, 8 for the corrections to core classes). They live in `TKAlgoExt/GTests/BRep*Ext_Test.cxx` (one file per package) with a shared
session and oracles in `BRepExtTest_Support.*`, in the `GeomAlgoExt` and `BRepAlgoExt` GTests, and in the test
files of the corrected classes `CSLib_Class2d`, `BRepTopAdaptor_FClass2d`, `BndLib`, `gp_Dir`, `gp_Dir2d` and
`BRepTools`. The Patch tests hold every fill they make to the same net at rest, turned, mirrored and at a
thousandth and a thousand times the size. For example:

```sh
OpenCascadeGTest --gtest_filter='BRep*Ext*:CSLib_Class2d*:BRepTopAdaptor*:BndLib*'
```

Every expected value comes from a closed form or a definition that does not use the code under test:
Steiner's volume of a box rounded on every edge; the rolling-ball definition of a blend (distance of the ball
centre to the input minus the radius); Archimedes' volume of two crossed cylinders; the frustum and Pappus
volumes of skins and sweeps; the circular segment of three arcs closed by a chord; the position of the circle
touching a line and an arc; and, for bridge curves, the curvature, change of curvature and torsion of a circle
and a helix at the ends. Each test names the defect it would catch.

Operations without a test of their own yet: deformation of whole bodies (`deformBody`, `deformSupports`),
Boolean imprints along faces' normals, and the curvature-continuous fills between plates meeting at three or
more edges (not yet supported). Of the tests that exist, the ones for regions, thickening, bridging and
aligning were each shown to fail with the code or the continuity order broken; the others have not been
mutation-checked.

## Repository layout

| Path | Content |
| --- | --- |
| `src/ModelingAlgorithms/TKAlgoExt/` | The nine operation packages and their GTests |
| `src/ModelingAlgorithms/TKGeomAlgo/GeomAlgoExt/` | Numerical utilities |
| `src/ModelingAlgorithms/TKMesh/BRepAlgoExt/` | Charts, trim integration, properties |
| `adm/cmake/occt_extensions.cmake` | Build options of the extension packages |
| `adm/scripts/extensions/` | Benchmarks, their results, the chart script, parity harness |
| `dox/extensions/` | Documentation of the additions and the images of this file |

## Upstreaming

OCCT would reasonably take this in pieces: (1) the corrections that apply to `master` today, each alone;
(2) the interval and Bernstein utilities; (3) the spline and rim fits with their GTests; (4) `BRepToolsExt` and
`BRepCurveExt`, which are self-contained and have the clearest closed-form checks; (5) `BRepOffsetExt`; (6) the
blend and Patch solvers last. The namespace `BRepExt`, the handle interface and the diagnostic hooks are
artefacts of how the code was first used. Moving the packages to OCCT's class conventions with
`Standard_Failure` is mechanical but changes behaviour (OCCT code catches `Standard_Failure`), so it is left for
review.

## Open work

Known limits and what remains:

1. **Blend solver.** A final recovery construction partitions joined contact strips at isolated
   shared-edge events, with matching pcurve restrictions and retained edge intervals; it keeps the
   fitted sheet unchanged and follows every previous construction. Overlapping or simultaneous
   boundary events still need a more general contact arrangement (see [operations.md](dox/extensions/operations.md#contact-ownership-on-a-joined-blend-sheet)).
   Some configurations where several blends meet are refused although a valid result exists:
   trimming a blend group needs an owner for every retained region and contact before chain ends are cut, and
   the centre path must be parametrised by arc length along the marched curve, not by the edge's parameter.
   Curvature-continuous and conic sections at corners where four or more blends meet are also refused. An edge whose faces meet at a
   steady angle is blended at any opening from 5 to 175 degrees, convex or concave (measured on planes,
   plane and cylinder, two cylinders). Where the angle turns along the edge so far that the balls' centre
   path bends more tightly than the ball (the rim of a hole leaning a with r >= (R + r) cos a: from about
   70 degrees), the rolling-ball surface crosses itself and is refused by name until it is trimmed there;
   short of that the rim is the balls' envelope. A final circuit route attaches a fixed-contact conic
   replacement on a cylindrical cavity and transverse plane, retaining the envelope outside its
   crossing sector. Its validity, closed-solid, sampled contact and retained-section checks match
   the other delivered blend routes; continuous spline and embedding certificates remain open
   (see the conic-window section in [operations.md](dox/extensions/operations.md)). Historical stored results over pairs of supports include
   values from an earlier construction; they do not establish that this solver delivered those results.
   Several curvature-continuous results recorded before the internal tangent-continuity check are also
   refused by that check. Recovery requires meeting the current seam and section-law bounds.
2. **Section law for curvature continuity.** `BRepBlendExt_Policy` supplies the conic section only. A published
   bound for a curvature-continuous conic or quintic section as a function of the support's curvature is
   missing.
3. **Fill.** The interior acceptance gate and the flat-support curvature limit have no derivation; across two
   rims of an open boundary the gate's reference is the Hermite bridge between them, which is reasoned, but
   the factor it allows over it is the same judgement, and gentle bulges where an elbow turns pass more often
   than against the rims alone. The ladder of nets stops where two nets in a row agree; its rungs are not
   nested spline spaces, so on some loops the distance between nets does not shrink (a plateau, or nets
   that never settle) and the stop is then a judgement, and a guide's miss can rise from one rung to the
   next (a side fill whose guide ends near a G2 edge is refused so). At a corner whose
   faces bend differently the bend asked eases to the one a surface can have over one span of the coarsest
   net, a judgement: G2 is met on the middle of each side only, the corner's miss is reported, and the coarse
   nets of such loops can come out lopsided. Exact G2 at a regular corner whose supports have
   different second fundamental forms is impossible: continuity would require the fill's one
   form there to equal both (the [corner compatibility proof](dox/extensions/operations.md#brepatchext-boundary-first-fill)).
   This includes a bending strip meeting a flat strip and rounds with unequal forms; local
   refinement or several internally G2 faces cannot remove the obstruction. The eased fill with
   its reported corner miss relaxes that impossible request. Some non-four-sided loops remain
   refused because the finite net misses the flat-support limit beyond the ease; the proof does
   not establish that this relaxed problem is impossible. Fills through
   three profiles with curvature continuity can miss the curvature-share limit. Open boundaries are chained by
   `BRepPatchExt_Chain` (nearest-end growth, gaps closed with straight segments); a general theory of which
   chain a user means is not attempted. A loop that is a height over no plane is filled over a faired
   disc laid flat conformally (Boundary First Flattening); the net over that chart is a tensor net, so a
   loop that needs fine spans in one place needs them right across. Its 640-pole forecast is a
   computation policy based on osculating-circle derivatives, not a general approximation bound or
   a necessary pole count. Complete curvature enforcement, converged geometric feasibility and
   support-derived collar initialization remain open. Not yet supported: curvature-continuous fills between plates
   meeting at three or more edges at once, strips with bending flows, flow-direction modes on steep
   supports, and profile-curve panels.
4. **Offsets and thickening.** Strict local focal crossings on non-periodic rational spline sheets are
   refused by interval witnesses; periodic supports, inconclusive bounds and global collisions remain
   unproven. A complete circular tube's outward parallel body closes its hole by revolving its clipped
   meridian. Offsetting solids whose blends meet the wall thickness exactly has no exactness proof. In the face offset
   the largest angle between the faces along a tangent edge is found by sampling, not bounded.
   Continuing an approximately tangent contact past its trim remains open. Endpoint gaps are checked
   independently of result tolerances, including both seam pcurves; chart cuts follow their
   junctions and closed nonperiodic section curves use independent endpoint stations.
   See [the endpoint certificate limits](dox/extensions/operations.md#face-offset-endpoint-certificate-limits).
5. **Speed.** Measured against stock followed by the same check: a single circular edge 1.05 times, the round
   blended into earlier rounds 1.19 times, the curve offset 1.5 times, the hollow 2.6 to 3.2 times, the smooth skin
   through square, circle and square 4.3 times. Everything else in the [table](dox/extensions/benchmarks.md) is level
   or faster. The test that a fitted blend does not cross itself is evaluated on every fitted blend; a curvature bound
   that proves the sections cannot run back would let it be skipped (not done). A profile of that skin
   ([benchmarks.md](dox/extensions/benchmarks.md#smooth-skin-through-square-circle-square)) puts a third of
   it in arc-length placement of points while a grid is refined, a quarter in sewing, a sixth each in face
   construction and in orienting the solid.
6. **Tests.** The list above.
7. **Conventions.** Formatting is OCCT's (`clang-format` 18.1.8 with the repository's file), headers carry the
   repository's licence text, comments use `//!`, and files are named `Package_Class`. Not done: the packages
   use free functions in namespace `BRepExt`, the handle interface and `std::runtime_error` failures where OCCT
   style would be classes with `IsDone` / `Shape`, `Standard_Failure` and `Message_ProgressRange`. Public
   free functions carry `Standard_EXPORT`; a Windows DLL build has not been tried.
8. **Policy.** The choices listed under [Section policy](#section-policy) are not mathematics and are not
   included.

## Licence and contribution

New files are Copyright (c) 2026 Timon Gentzsch under OCCT's terms: GNU LGPL version 2.1 with the OCCT
exception (`LICENSE_LGPL_21.txt`, `OCCT_LGPL_EXCEPTION.txt`). Existing OCCT files keep their notices.

OCCT's own [contribution workflow](dox/contribution/contribution_workflow/contribution_workflow.md) asks
contributors to sign its Contributor License Agreement before a change is accepted upstream, and to title a
change `<Group> - <Summary>`. The commits here follow that title form. No CLA has been signed; it would be
needed before offering any of this to OCCT.

---

# Open CASCADE Technology

Open CASCADE Technology (OCCT) is a software development platform providing services for 3D surface and solid modeling, CAD data exchange, and visualization. Most of OCCT functionality is available in the form of C++ libraries. OCCT is ideal for developing software dealing with 3D modeling (CAD), manufacturing/measuring (CAM), or numerical simulation (CAE).

## License

Open CASCADE Technology is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License version 2.1 as published by the Free Software Foundation, with a special exception defined in the file `OCCT_LGPL_EXCEPTION.txt`. Consult the file `LICENSE_LGPL_21.txt` included in the OCCT distribution for the complete text of the license.

Alternatively, Open CASCADE Technology may be used under the terms of the Open CASCADE commercial license or a contractual agreement.

**Note:** Open CASCADE Technology is provided on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND. The entire risk related to any use of the OCCT code and materials is on you. See the license text for a formal disclaimer.

## Packaging

You can receive certified versions of OCCT code in different packages:

- **Snapshot of Git repository:** Contains C++ header and source files of OCCT, documentation sources, build scripts, and CMake project files.
- **Complete source archive:** Contains all sources of OCCT, generated HTML and PDF documentation, and ready-to-use projects for building on all officially supported platforms.
- **Binary package (platform-specific):** In addition to the complete source archive, it includes binaries of OCCT and third-party libraries built on one platform. This package allows using OCCT immediately after installation.

Certified versions of OCCT can be downloaded from:
- [Open CASCADE Releases](https://dev.opencascade.org/release)
- [GitHub Releases](https://github.com/Open-Cascade-SAS/OCCT/releases)

You can also find OCCT pre-installed on your system or install it from packages provided by a third party. Note that packaging and functionality of such versions can be different from certified releases. Please consult the documentation accompanying your version for details.

## Documentation

Documentation is available at the following links:
- [Latest version](https://dev.opencascade.org/doc/overview)
- [Version 8.0](https://dev.opencascade.org/doc/occt-8.0.0/overview)

Documentation can be part of the package. To preview documentation as part of the package, open the file `doc/html/index.html` to browse HTML documentation.

If HTML documentation is not available in your package, you can:

- **Generate it from sources:** You need to have Tcl and Doxygen 1.8.4 (or above) installed on your system and accessible in your environment (check the environment variable PATH). Use the batch file `adm/gendoc.bat` on Windows or the Bash script `adm/gendoc` on Linux or OS X to (re)generate documentation.
- **Generate together with sources:** You need to have CMake and Doxygen 1.8.4 (or above) installed on your system. Enable `BUILD_DOC_Overview` CMake parameter and set the path to Doxygen `3RDPARTY_DOXYGEN_EXECUTABLE`. Then build ALL or only `Overview`.
- **Read documentation in source plain text (Markdown) format** found in the subfolder `dox` or [GitHub Wiki](https://github.com/Open-Cascade-SAS/OCCT/wiki).

See [dox/build/build_documentation/building_documentation.md](dox/build/build_documentation/building_documentation.md) or [Building Documentation](https://dev.opencascade.org/doc/occt-8.0.0/overview/html/build_upgrade__building_documentation.html) for details.

## Building

In most cases, you need to rebuild OCCT on your platform (OS, compiler) before using it in your project to ensure binary compatibility.

Consult the file [dox/build/build_occt/building_occt.md](dox/build/build_occt/building_occt.md) or [Building OCCT](https://dev.opencascade.org/doc/overview/html/build_upgrade__building_occt.html) or [Building OCCT Wiki](https://github.com/Open-Cascade-SAS/OCCT/wiki/build_upgrade) for instructions on building OCCT from sources on supported platforms.

## Version

The current version of OCCT can be found in the file [`adm/cmake/version.cmake`](adm/cmake/version.cmake).

## Development

### Bug Tracker
- [GitHub Issues](https://github.com/Open-Cascade-SAS/OCCT/issues)
- [OCCT Tracker](https://tracker.dev.opencascade.org/)

For information regarding OCCT code development, please consult the official OCCT Collaborative Development Portal:
- [OCCT Development Portal](https://dev.opencascade.org)

### Forum and Discussions
- [OCCT Forums](https://dev.opencascade.org/forums)
- [GitHub Discussions](https://github.com/Open-Cascade-SAS/OCCT/discussions)

