# Corrections to existing OCCT algorithms

Twenty-one changes to existing OCCT code, one commit each, on top of OCCT `b8f597c6` (30 July 2026).
Each entry states the problem, the derivation that makes the change correct, and its limits.
Derivations are elementary arguments, not machine-checked proofs, unless a paper is cited.
Existing OCCT tolerances and limits are kept unless an entry names an added constant.

To review one change alone: `git switch -c review b8f597c6 && git cherry-pick COMMIT`.

"Upstream" gives the state when each was compared with OCCT `master` at `53706f78a8`:
*applies* means it cherry-picks onto `master` and the toolkit builds; *changed upstream* means
`master` has since redesigned that code, so the change is not ported but the problem statement
stands. All are independent, except that offset tolerance cancellation needs the offset progress
guard (`Geom_OffsetSurfaceUtils-cancellable-progress`) first.

| Entry | Commit | Kind | Upstream |
| --- | --- | --- | --- |
| [Approx_BSplComputeLine-local-endpoint-parameters](#approx_bsplcomputeline-local-endpoint-parameters) | `806f921d8a` | wrong parameter scale | applies |
| [BOPAlgo_RemoveFeatures-extend-across-feature](#bopalgo_removefeatures-extend-across-feature) | `af44c0a702` | needless work | applies |
| [BRepAlgo_Loop-sort-vertices](#brepalgo_loop-sort-vertices) | `2c9b2590d8` | quadratic loop | applies |
| [BRepCheck_Solid-stop-after-enclosed-region](#brepcheck_solid-stop-after-enclosed-region) | `561fa613b4` | uninitialised read | applies |
| [BRepFill_Draft-orientation-last-sample](#brepfill_draft-orientation-last-sample) | `2401e8b921` | out-of-range query | applies |
| [BRepGProp_Gauss-used-errors](#brepgprop_gauss-used-errors) | `bcf74fe6c5` | unbounded scan | changed upstream |
| [BRepMesh_Delaun-hilbert-insertion](#brepmesh_delaun-hilbert-insertion) | `3fa3feb2d4` | insertion order | applies |
| [BRepMesh_DelaunayDeflectionControlMeshAlgo-node-normals](#brepmesh_delaunaydeflectioncontrolmeshalgo-node-normals) | `e80dd7d596` | repeated work | changed upstream |
| [Geom_OffsetSurfaceUtils-cancellable-progress](#geom_offsetsurfaceutils-cancellable-progress) | `6593918e73` | missing cancellation | applies |
| [BRepOffset_MakeOffset-cancellable-tolerances](#brepoffset_makeoffset-cancellable-tolerances) | `bdd6b9a822` | missing cancellation | applies (needs the guard) |
| [BRepOffset_MakeOffset-input-order](#brepoffset_makeoffset-input-order) | `b045680776` | address-dependent order | applies |
| [BRepTopAdaptor_FClass2d-wire-bounds](#breptopadaptor_fclass2d-wire-bounds) | `c8c90a2c33` | needless work | changed upstream |
| [Extrema_ExtPExtS-spline-across](#extrema_extpexts-spline-across) | `a71c466051` | missed separability | applies |
| [IntTools_BeanFaceIntersector-whole-range](#inttools_beanfaceintersector-whole-range) | `dc671c167b` | missed intersections | applies |
| [IntTools_EdgeFace-line-plane-touch](#inttools_edgeface-line-plane-touch) | `def5b1740c` | missed intersections | applies |
| [IntTools_FaceFace-plane-across-sweep](#inttools_faceface-plane-across-sweep) | `5ce72fca97` | lost exact construction | applies |
| [TopOpeBRepBuild_WireEdgeClassifier-box-apart](#topopebrepbuild_wireedgeclassifier-box-apart) | `14f92b09b4` | needless work | applies |
| [gp_Dir-restored](#gp_dir-restored) | `ac6046d3d7` | lost bits on reload | not assessed |
| [BndLib-torus-patch-box](#bndlib-torus-patch-box) | `810a8e12a5` | box misses the patch | not assessed |
| [CSLib_Class2d-v-bucket-index](#cslib_class2d-v-bucket-index) | `341272613f` | linear scan per query | not assessed |
| [BRepTopAdaptor_FClass2d-periods-read-once](#breptopadaptor_fclass2d-periods-read-once) | `6e6ea43f94` | repeated work | not assessed |

A twenty-second entry, [BRepSkinExt-section-vertices-exact](#brepskinext-section-vertices-exact), fixes a defect
in the extension package `BRepSkinExt` rather than in OCCT.

### gp_Dir-restored

Commit `ac6046d3d7`. Upstream: not assessed.

Problem: the BREP readers (`BinTools_CurveSet`, `BinTools_Curve2dSet`, `BinTools_SurfaceSet` and the
`GeomTools` sets, binary and text) rebuilt every direction with `gp_Dir::SetCoord`, which keeps the coordinates
only when their computed square length is 1 within `gp::Resolution()` (about 1e-308), so in effect only for exact
unit vectors, and divides by the norm otherwise. A direction normalized once has a computed square length a few
units of roundoff off 1, so a body read back from a file differed from the one written in the last bits of its
oblique lines, line pcurves and extrusion directions. Measured: of 1e5 random directions, 32 per cent (3D) and 31
per cent (2D) changed when made again from their own coordinates; a body and the same body saved and loaded
meshed with the same nodes and different triangles (an incremental mesh comparison, seed 27).

Derivation and limits: for v normalized as y_i = fl(v_i / D), D = fl(sqrt(fl(sum v_i^2))), with unit roundoff
u = 2^-53 and g_n = n u / (1 - n u): fl(sum v_i^2) = |v|^2 (1 + t1), |t1| <= g3; D = |v| (1 + t1)^(1/2) (1 + e0),
|e0| <= u; y_i = (v_i / D)(1 + e_i), |e_i| <= u. So |y|^2 lies within 2u + g3 + 2u of 1, and the square length
computed from y, |y|^2 (1 + t2) with |t2| <= g3, within 10 u + O(u^2) of 1 (8 u in the plane); contracting
products into FMA only removes roundings. `gp_Dir::Restored` and `gp_Dir2d::Restored` keep the coordinates where
that computed square length is within `gp::UnitSquareBand()` = 8 DBL_EPSILON = 16 u of 1, and normalize as before
otherwise, and the six readers use them. Every direction `gp_Dir` or `gp_Dir2d` normalized is therefore read back
bit for bit, and restoring is idempotent. Normalization itself is unchanged: widening `gp_Dir`'s own band instead
changed last bits throughout the kernel and made 32 of 90 runs of the incremental mesh fuzz fail. The band is the
one added constant; 16 u against the 10 u bound leaves the margin for directions made by other arithmetic within a
few units of roundoff of unit length. Not covered: an axis placement's X and Y directions, which the readers rebuild
from its main and X directions by cross products (`gp_Ax3(P, N, Vx)`), and come back exact only where those
products are (operations.md, open). Tests: random vectors of sizes 2^-400 to 2^400 normalize within 10.5 u (8.5 u)
of unit square length, are restored bit for bit, and vectors 24 units of roundoff off unit length restore
idempotently; 200 random oblique lines, plane normals and lines on the planes read back bit for bit through the
binary and the text format (the test fails with the previous readers).

### BndLib-torus-patch-box

Commit `810a8e12a5`. Upstream: not assessed.

Problem: `BndLib::Add(gp_Torus, UMin, UMax, VMin, VMax, ...)` returns a box that does not contain the patch. It boxes the arcs
that points of the section circle, every 45 degrees, trace about the axis; those points lie on the circle, so the part of the
circle between two of them lies outside the polygon through them. Random tilted tori: 0.6 per cent of 7.4 million surface
points fall outside, by up to 7.1 per cent of the minor radius (the inscribed octagon's bound is 1 - cos(pi/8) = 7.6 per cent).

Derivation and limits: the patch is the rotation about the axis, through the U range, of an arc of the section circle. That arc
lies in the convex hull of its end points and, for each sub-arc of at most 45 degrees, of the intersection of the tangents at
the sub-arc's ends, which is on the bisector at Ri / cos(step / 2) from the section's centre. A linear functional on the
rotation of a convex set about an axis is largest at a vertex of the set, so the union of the arcs' boxes that these points
trace contains the patch. Rotation circles are boxed as before (by the circumscribed octagon, so a full circle of radius R is
boxed to R / cos(pi/8)); the degenerate torus (R < r) keeps its own routine. Tests: sampled surface points of tilted tori
with reversed and long V ranges lie in the box (the test fails on the previous code), and the inner half of a torus does
not reach the outer equator.

### CSLib_Class2d-v-bucket-index

Commit `341272613f`. Upstream: not assessed.

Problem: `CSLib_Class2d::SiDans` scans every vertex of the polygon for each query, and up to five times: once for the point
and once for each corner of its tolerance box. `BRepTopAdaptor_FClass2d::Perform` calls it for every wire of a face, so
classifying points against a face with a rounded boundary costs the number of boundary samples each time. A profile of a blend
into earlier rounds spent a quarter of the call in it (68 000 queries, about 4 000 instructions each).

Derivation: what an edge can contribute to a classification at height Y is a crossing of the ray (needs Y between its end
heights), a vertex within the tolerance box (needs Y within the V tolerance of the vertex's height) or the interpolated edge
passing within the V tolerance (the interpolated height lies between the end heights, so Y is within the tolerance of the
edge's height range). So an edge matters only when Y lies in its V range widened by the tolerance. The edges of a polygon of
32 vertices or more are listed once under every bucket (of about four edges each) that their widened range meets, and a query
reads the one bucket of its Y. The same edges are examined in the same order with the same arithmetic, so every answer is the one
the scan gives. Test: random wavy polygons of 12 to 4 096 vertices and a staircase whose vertices lie on bucket boundaries,
queried at random points, at vertices, a fraction of the tolerance from vertices and at edge middles; the three answers (inside,
outside, uncertain) are compared with a reference written from the definition, which itself is checked against the unchanged
class. Effect: a blend into earlier rounds 7.96 to 7.16 ms, two crossed cylinders 12.9 / 15.4 to 8.4 / 10.3 ms.

### BRepTopAdaptor_FClass2d-periods-read-once

Commit `6e6ea43f94`. Upstream: not assessed.

Problem: `Perform` and `TestOnRestriction` build a `BRepAdaptor_Surface` of the face on every call to ask whether the surface
is periodic in U and V and by how much. The answer cannot change after construction (about 6 per cent of a blend into earlier
rounds). The four numbers are now read in the constructor. Test: points a whole number of periods from a half cylinder's wall
classify as the wall does, in and out, before and after the change.

### BRepSkinExt-section-vertices-exact

Commit `6efcd84dcc` (extension package, not OCCT code). Found by the skin rows of the [benchmark](benchmarks.md).

Problem: a smooth skin through a square, a circle and a square was an unconnected solid (`BRepCheck_NotConnected`, three
pairs of coincident unsewn edges at the square's corners) where stock `BRepOffsetAPI_ThruSections` made a valid one. The
circle's arcs are rational and the square's edges are not, so the arcs are made again without weights to put every edge of every
section on one basis; neighbouring arcs then ended 1.4e-6 apart, above the sewing tolerance of 1e-6, and the corner curves of the
faces on either side did not coincide.

Fix and limit: the end poles of a clamped curve are its end points, which are exact on the curve being replaced, so the
replacement's end poles are set to them. The fit moves by at most the error it already had (1.4e-6 at the ends), and the
sections stay on the result to 7e-15. Test: the skin is a valid solid and every sampled point of the three sections lies on its
faces within the sewing tolerance (the test fails without the change).

### Approx_BSplComputeLine-local-endpoint-parameters

Commit `806f921d8a`. Upstream: applies; native toolkit build passed.

Problem: Endpoint fitting previously mixed physical parameters with normalized local coordinates.

Derivation and limits: For a positive local span h=b-a, s=(t-a)/h and dC/dt=(dC/ds)/h. Solve the endpoint Bezier fit on [0,1] and rescale both endpoint tangents by h. Assumes strictly increasing endpoint parameters. No new threshold or sampling count.

### BOPAlgo_RemoveFeatures-extend-across-feature

Commit `af44c0a702`. Upstream: applies; native toolkit build passed.

Problem: Extending both directions unnecessarily enlarges adjacent face domains during feature removal.

Derivation and limits: On a forward face an oriented parameter tangent (du,dv) has outward transverse direction (dv,-du). Extend only coordinates with a nonzero outward component; a missing pcurve retains the full extension. This local derivation does not prove that five samples capture every tangent on a curved edge. Four subintervals are an unproven sampling choice. Existing angular and geometric resolutions are retained.

### BRepAlgo_Loop-sort-vertices

Commit `2c9b2590d8`. Upstream: applies; native toolkit build passed.

Problem: Repeated parameter evaluation and adjacent exchange sorting take quadratic work.

Derivation and limits: Deduplicate in encounter order, evaluate each remaining parameter once, then stable-sort it. For finite parameters this preserves the stable order of the adjacent exchange sort, including equal parameters, while reducing sorting comparisons from quadratic to O(n log n). Assumes the shape map equality agrees with the old duplicate comparison. No new geometric constants.

### BRepCheck_Solid-stop-after-enclosed-region

Commit `561fa613b4`. Upstream: applies; native toolkit build passed.

Problem: Continuing after an enclosed-region error reaches classifiers that were never initialized.

Derivation and limits: Once a second growth shell adds EnclosedRegion, later classifiers have not been initialized. Return before their use. The already-recorded invalidity is retained; no claim is made about collecting further diagnostic errors. No new constants.

### BRepFill_Draft-orientation-last-sample

Commit `2401e8b921`. Upstream: applies; native toolkit build passed.

Problem: Rounded final abscissae can exceed the wire length and query an absent location law.

Derivation and limits: Exact sample positions satisfy 0 <= (i-1)*length/N <= length. Floating-point multiplication can place the endpoint outside this closed interval. min(sample,length) restores the endpoint bound and prevents querying an absent location law. Interior arithmetic and sample count are unchanged. No new constants.

### BRepGProp_Gauss-used-errors

Commit `bcf74fe6c5`. Upstream: not ported: changed upstream implementation.

Problem: Adaptive integration searched unwritten zero-initialized error capacity.

Derivation and limits: The error vector has an explicitly zero-initialized suffix. Its maximum is the prefix maximum unless that is negative, in which case the first suffix zero wins. Scan the used prefix with the same strict comparison and reproduce that suffix case. This retains the selected integration interval and arithmetic. No new numerical constants.

### BRepMesh_Delaun-hilbert-insertion

Commit `3fa3feb2d4`. Upstream: applies; native toolkit build passed.

Problem: Geometrically scattered insertion order makes mesh refinement expensive.

Derivation and limits: Map the parameter box to a 16-bit square grid and order nodes by recursive quadrant rotation, breaking ties by node index. This is an ordering change only; deflection controls are unchanged. Locality is a performance heuristic, not a runtime theorem; the idea follows [Skilling 2004](references.md). The 16-bit grid (65535 maximum, loop starts at bit 15) is an unproven resolution choice; triangle connectivity and sampling may change.

### BRepMesh_DelaunayDeflectionControlMeshAlgo-node-normals

Commit `e80dd7d596`. Upstream: not ported: changed upstream implementation.

Problem: Repeated normal estimation at unchanged mesh nodes duplicates work.

Derivation and limits: For an unchanged surface and identical parameter pair, NormEstim is the same deterministic function call. Cache its value and defined flag keyed by node index plus exact parameter coordinates, invalidating between optimizations. Coordinate validation prevents stale results after index reuse. No tolerance, angular limit or normal formula changes.

### Geom_OffsetSurfaceUtils-cancellable-progress

Commit `6593918e73`. Upstream: applies; native toolkit build passed.

Problem: Offset normal retries can stop making representable progress without observing cancellation.

Derivation and limits: An iteration with unchanged stored parameters evaluates the same normal again and cannot advance. Reject nonfinite or nonrepresentable steps without an iteration bound; retain arithmetic for progressing steps. A borrowed thread-local Message_ProgressRange is scoped by RAII, so nested evaluations restore their callers even on exceptions. No new numerical constants.

### BRepOffset_MakeOffset-cancellable-tolerances

Commit `bdd6b9a822`. Upstream: applies textually; requires progress-guard dependency to build.

Problem: Edge tolerance sampling lacked cancellation and loaded derivative-dependent adaptors unnecessarily.

Derivation and limits: Positional D0 evaluation uses only the loaded pcurve and surface, not the endpoint derivative adaptors. Load curve then surface and retain all positional samples and arithmetic. Check the progress range before applying incomplete tolerance results. Depends on the offset evaluation progress guard. The new progress weight 1 is an unproven presentation share, not a geometry bound or deadline.

### BRepOffset_MakeOffset-input-order

Commit `b045680776`. Upstream: applies; native toolkit build passed.

Problem: Offset root traversal depended on allocation-sensitive map order.

Derivation and limits: Sort map entries by TopExp input index before building their generated faces. The traversal becomes a function of input topology instead of allocation addresses, provided every key has an input index. Stable sorting retains existing order for equal indices. No tolerance or geometric constant changes.

### BRepTopAdaptor_FClass2d-wire-bounds

Commit `c8c90a2c33`. Upstream: not ported: changed upstream implementation.

Problem: Using a face-wide box for each wire prevents rejecting unrelated holes cheaply.

Derivation and limits: A polygon lies in the coordinate-wise extrema of its vertices. A query outside that box beyond the classifier tolerance cannot be inside the polygon. Use that wire box for quick rejection while retaining the identical polygon ray test inside and the old fallback for zero-area bounds. No new constants.

### Extrema_ExtPExtS-spline-across

Commit `a71c466051`. Upstream: applies; native toolkit build passed.

Problem: A generic surface extremum solve ignores separability of transverse spline extrusions.

Derivation and limits: For unit D and C(u).D=h constant, S(u,v)=C(u)+vD gives squared distance ||P-h_proj-C(u)||^2+(v-v0)^2. The interior minimum therefore separates into a curve projection and v=v0. Keep the generic path for boundary or unresolved cases. Existing capacity of four extrema is retained: this does not promise every boundary or non-nearest extremum. Existing tolerances are unchanged.

### IntTools_BeanFaceIntersector-whole-range

Commit `dc671c167b`. Upstream: applies; native toolkit build passed.

Problem: Localized edge/face searches may omit intersections outside their selected neighborhood.

Derivation and limits: Use the existing whole-range search for all surface kinds instead of choosing a localized search for certain polynomial surfaces. Both seek the same zero distance locus; no theorem guarantees completeness of the numerical global search. This preserves the previously selected implementation, not a new intersection algorithm. No new constants.

### IntTools_EdgeFace-line-plane-touch

Commit `def5b1740c`. Upstream: applies; native toolkit build passed.

Problem: Numerical touch detection can miss the analytic interior intersection of a line and a plane.

Derivation and limits: For L(t)=p+td and plane (x-q).n=0, t=(q-p).n/(d.n) when d.n is nonzero. Use that intersection only inside both trimmed parameter ranges, retaining the generic extrema search elsewhere. Angular nonparallel recognition uses the existing Precision value. No new thresholds.

### IntTools_FaceFace-plane-across-sweep

Commit `5ce72fca97`. Upstream: applies; native toolkit build passed.

Problem: The generic section search loses the exact translated-curve construction of transverse extrusions.

Derivation and limits: For S(u,v)=C(u)+vD, D parallel to plane normal n, and C(u).D constant, v=(q-C(u0)).n/(D.n) is constant. Translate the basis curve by vD; its extrusion pcurve is (u,v). Preserve boundedness for section endpoints, and retain the general path near trimmed ends. Assumes positive-weight spline poles for the constant-height hull argument. Existing OCCT tolerances only.

### TopOpeBRepBuild_WireEdgeClassifier-box-apart

Commit `14f92b09b4`. Upstream: applies; native toolkit build passed.

Problem: Wire classification performs expensive curve tests even for disjoint parameter boxes.

Derivation and limits: A loop is contained in the bounds of all its existing pcurves. If two expanded boxes are disjoint, one loop is in the unbounded complement of the other box and can reuse its outside classification. Missing pcurves and block classification keep the old path. The outside probe uses xmax+width+1; 1 is an unproven positive displacement, not an error bound. UV expansion retains stored edge tolerances and therefore inherits their unit assumptions.
