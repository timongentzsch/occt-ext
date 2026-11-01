<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# OCCT mathematical extensions and corrections

This local review branch starts at official OCCT `b8f597c6`, followed by
`055a9a8a` (WebAssembly exception conversion) and `c7f60ffe` (null curve guard).
It preserves OCCT history. The 17 corrections below are separate commits;
new finite mathematical utilities form package `GeomAlgoExt` in toolkit `TKGeomAlgo`.
This is a partial separation of application algorithms, not an assertion that
all application solvers have become independent OCCT packages.

## Exact-section blend fitting

`GeomAlgoExt_BlendFit` adds a tensor-surface fit to caller-supplied exact rational
sections and support contacts. A degree-p B-spline is interpolated along the rail
at Greville sites, or supplied terminal Hermite jets, for each homogeneous pole
coordinate and both support parameter curves. For X=wP, differentiation gives
P'= (X'-Pw')/w and P''=(X''-2P'w'-Pw'')/w. For a support contact (U,V),
the second fundamental form applied to its velocity is
n dot (S_uu U'^2+2S_uv U'V'+S_vv V'^2); dividing by squared contact speed
compares normal curvatures. Basis and derivative identities use de Boor [3].

Assumptions: finite valid clamped bases, consistent section pole count, positive
weights, callbacks returning the documented row width 4S+4 and valid derivatives,
unit support normals, and live input/callback lifetimes. Collocation can be
singular. Sampled inter-site positions and curvature are **not global error or
G2 certificates**. Stability, the stall/refinement rules and flat-support limits
are not proved. No section construction law or measured coefficient is supplied
by this module; values and section/regime/contact callbacks are explicit inputs.

`BlendSectionFitter::Advance` owns one original numerical refinement step.
A caller may retain its existing outer polling loop, without moving a host callback
into the numerical code. The convenience `fitBlendSections` accepts OCCT progress
for native clients. This separation preserves the old app polling path: an opaque
range check did not meet the zero-increase hot polling gate. It does not claim that
the rest of the rail marching, section construction, topology or attachment moved.

Retained **unproven choices**: maximum 16 refinements, 4000 longitudinal poles,
curvature refinement only before round 8, three interior probes at 1/4,1/2,3/4,
late-break discovery at rounds 2 and 3, 128 regime scan intervals, 60 bisections,
three-round stall history, factor 1000 and halving criterion. Degree, bend budget,
merge/degeneracy/length tolerances and law decisions are supplied explicitly.
The quarter-cylinder GTest checks the defining circle, axial coordinate and
contact curves at three scales; a separate test checks cancellation before sampling.

`adm/scripts/extensions/module-benchmark.sh blend-fit FORK_BUILD STOCK_BASE_BUILD`
compares the same quarter-cylinder on the vertical (0,0,z) edge of a 10s by 10s by
20s box, radius s, s=0.1,1,10. The new method fits the exact rational sections;
stock `BRepFilletAPI_MakeFillet` builds the entire filleted body. Setup/edge selection
and a 33x33 defining-circle/axial survey are outside the clock. This is deliberately
unequal topology work, **not an end-to-end fillet speedup**.

| Scale | Stock whole-body CPU ms | Section-fit CPU ms | Stock circle distance | Fit circle/axial distance |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 0.3025 | 0.0120 | 2.776e-17 | 8.882e-16 |
| 1 | 0.2950 | 0.0110 | 4.441e-16 | 1.066e-14 |
| 10 | 0.2955 | 0.0110 | 3.553e-15 | 8.527e-14 |

Ten successes/zero refusals after one warmup per row; process CPU, GCC12.2 -O3,
TBB off, Ryzen 9 7950X3D shared host. Stock is more accurate on these finite
fixtures. Min/max values are in [blend-fit-measurements.json](adm/scripts/extensions/blend-fit-measurements.json).
The independent native tests and timings were repeated after separating the numerical
step from the convenience loop. Integrated app CPU measurements separately gate
the behavior-preserving relocation.

## Boundary-first patch chart and interior survey

Package `BRepAlgoExt` uses existing toolkit `TKMesh`, because chart construction
requires its Delaunay mesher; it also uses existing lower topology/geometry tools.
`BRepAlgoExt_Chart` forms cotangent stiffness L, solves the interior equation
L_II u=K_I with zero boundary log scale, and obtains boundary turns from
pi-angleSum-(Lu)_boundary. For directions t_i and original lengths l_i,
minimize sum (delta l_i)^2/l_i subject to sum (l_i+delta l_i)t_i=0.
Lagrange multipliers give delta l_i=-l_i t_i^T G^-1 sum l_j t_j,
G=sum l_i t_i t_i^T, the implemented 2x2 closure. Harmonic x/y extension
then solves the Dirichlet energy equations. This follows Rohan Sawhney and
Keenan Crane (2017), [Boundary First Flattening](https://www.cs.cmu.edu/~kmcrane/Projects/BoundaryFirstFlattening/paper.pdf),
read sections 3–5, especially equations 6, 14 and 20. Their section 4.4 distinguishes
boundary interpolation from harmonic conjugation: this implementation uses
interpolation, so **exact discrete conformality is not guaranteed**.

Assumptions: finite regular surface, connected nondegenerate disc mesh, usable
boundary and nonsingular closure, positive resulting edge lengths and triangle
orientation. Cotangent Dirichlet systems are positive definite with fixed nonempty
boundary under the usual finite-element assumptions; arbitrary inputs can refuse.
The node cap, finite curvature survey and parameter-space lattice do **not** certify
3D deflection or global injectivity. The chord estimate k*l²/8 is local asymptotic,
not a proved global surface error bound. The chart may extrapolate off its mesh.

`BRepAlgoExt_Interior` keeps the original sampled curvature, fold and knot-step
surveys. Positive rational Bernstein weights put a boundary arc in its control
hull; outside the control-hull tubes, replacing arcs with chords preserves parity
under valid closed-wire/homotopy assumptions. Floating construction allowances
are unproved. Ambiguous trim points use OCCT classification. The sampled interior
maximum is **not** a curvature supremum. The caller supplies UV tolerance and the
bend ratio; the inherited ratio four is a host acceptance heuristic, with no native
default or theorem. The flat G2 limit is also still unproven and outside this package.

Retained **unproven choices**: 32000 chart nodes; CG residual target 1e-13 and
4n passes; minimum two edge subdivisions, half-spacing rim margin, and the
sqrt(triangleCount/2) lookup-grid sizing. Triangular-lattice sqrt(3)/2 and
sqrt(3)/4 are derived from equilateral geometry; cotangent 1/2 and triangle 1/2
are derived energy/area factors. Interior surveys use 40,64,100,300 samples;
Bezier tube goal 1e-4 span, squared-length cutoff 1e-24, recursion 12,
1e-10 span/length guards, 1e-12 length/collinearity/derivative cutoffs,
0.03 span exclusion and 1e-9 one-sided knot offset. These are inherited heuristics,
not tolerances inferred from the paper. No value or arithmetic order changed.

Native tests independently check planar pair distances and sphere curvature at
three scales, plus the explicit ratio gate. Run
`adm/scripts/extensions/module-benchmark.sh patch-chart FORK_BUILD STOCK_BASE_BUILD`
after building static Release TKMesh in both revisions. A rectangle and quarter
cylinder are the same surfaces/trims at scales 0.1,1,10. Stock provides meshing
and scaled carrier UVs; it has no equivalent boundary-first chart operation here.
That reference does less work and generates different meshes. The survey rows
compare stock GeomLProp's same 40x40 curvature grid with the trim-aware interior
survey. Geometry/angle checks are outside the process CPU clock.

| Scale 1 fixture | Stock CPU ms | New CPU ms | Stock/new max cosine error | Stock/new nodes |
| --- | ---: | ---: | --- | --- |
| plane-chart | 0.0760 | 0.2770 | 0 / 5.54e-13 | 4 / 75 |
| cylinder-chart | 0.1005 | 0.4510 | 0.000155 / 0.00222 | 18 / 101 |
| plane-curvature-survey | 0.0830 | 0.0965 | 0 / 0 | 0 / 0 |
| cylinder-curvature-survey | 0.1150 | 0.1555 | 0 / 0 | 0 / 0 |

Ten successes/zero refusals per row, after one warmup; GCC12.2 -O3, TBB off,
Ryzen 9 7950X3D shared machine. Trim-aware surveys are slower than the stock
untrimmed grid; the chart can also be slower than meshing/carrier UV alone.
This does not compare equivalent flattening solvers. All scales, extrema,
relative edge distortions and curvature errors are in
[patch-chart-measurements.json](adm/scripts/extensions/patch-chart-measurements.json).

## Constrained tensor-spline fitting and certified compaction

`GeomAlgoExt_SplineFit`, `GeomAlgoExt_RimFit` and `GeomAlgoExt_CompactSurface`
extend `TKGeomAlgo`; `GeomAlgoExt_Band` is their banded Cholesky engine.
They accept OCCT geometry and ordinary C++ containers. Refusals use
`GeomAlgoExt_Failure` with an unchanged stable reason/message; cancellation uses
`Message_ProgressRange`. An optional stage observer is a function pointer with
no clock, global state or application runtime. There is no cancellation added
inside the original fitting loops that had none.

For a tensor spline S(u,v)=sum P_ij N_i(u) M_j(v), scalar conditions are
linear forms a_r^T p=b_r in the free pole coordinates. The energy is
integral [alpha (|S_uu|²+2|S_uv|²+|S_vv|²) +
beta (|S_uuu|²+3|S_uuv|²+3|S_uvv|²+|S_vvv|²)] du dv.
The coefficients 1,2,1 and 1,3,3,1 count ordered mixed derivatives, hence
are derived. Substituting the basis gives a positive semidefinite quadratic
p^T H p+2g^T p+c. Weighted rows add A^T W A; eliminating fixed poles adds
their exact contribution to the right hand side. Cholesky solves the resulting
symmetric positive definite system when it is positive definite. Multipliers
correct row residuals using products with A M^-1 A^T. Stopped iterations are
retained exactly; convergence, consistency or global optimality are not promised.
The semidefinite claim assumes nonnegative energy coefficients, row weights and
density, finite regular knot spans and consistent constraint dimensions. Solvability
requires positive definiteness after fixed-pole elimination; these are caller
preconditions, not a proof that arbitrary inputs produce a usable surface.
Gauss quadrature with p+1 points integrates the polynomial basis products on
each knot span exactly, but arbitrary supplied density functions have no such
exactness guarantee. Basis locality/derivative identities follow de Boor (1972),
[On Calculating with B-Splines](https://ftp.cs.wisc.edu/Approx/bsplcalc.pdf),
DOI [10.1016/0021-9045(72)90080-9](https://doi.org/10.1016/0021-9045(72)90080-9),
read for this extraction. These algebraic derivations do not certify conditioning.

Rim correction linearizes the nearest-point displacement along a supplied unit
vector and solves min_x ||A x-r||²+d²||x||². Its normal equations are
(A^T A+d²I)x=A^T r; unobserved poles stay fixed. The nonlinear projection,
refinement decisions, summation order and refusal messages are preserved.
Compaction proposes OCCT knot removals, converts original and candidate to a
common piecewise Bernstein basis with outward intervals, and bounds displacement
by max_i ||Delta P_i||_2. Nonnegative basis functions summing to one imply this
bound over the whole rectangle, provided both conversions are enclosed, degrees
and full knot sequences agree, weights are nonrational, and the input is finite,
clamped and nonperiodic. No normal, trim or continuity bound follows from it.
Cancellation must escape the optional proposal's catch of OCCT builder failures.

All retained non-derived choices are explicit: fit diagonal holds 1e-13, factor
100 and cap 1e-6; stall count 12, stall share 0.02, maximum passes 400 and energy
ratio stop 1e-20; basis capacity degree 5 and five terms per condition; rim Newton
steps 8, determinant cutoff 1e-12, squared-step cutoff 1e-28, sliver ratio 0.2,
three solves, damping 0.03,
slow-progress ratio 0.5, cosine cutoff 0.17; compaction 512 trials, 16 removals,
262144 coefficients, 8*1024*1024 visits, proposal multiplier 1024, 33-point
rejection screen and 1e-10 screen allowance. These are implementation limits or
heuristics, not theorems, and none was tuned by this move. In particular 1e-20
is not binary64 epsilon squared. Cancellation polling cadences 1024 rows in rim
fitting and 16/256 items in compaction are inherited scheduling choices, not
proved latency bounds. The parent construction's interior acceptance
heuristic and flat G2 limit remain outside these modules and remain unproven.

The WASM SIMD path is retained. A portable two-lane path rounds each operation
through volatile binary64 stores to prevent FMA contraction. The finite parity
harness `adm/scripts/extensions/band-parity.cxx` exercises paired and ordinary
factors/solves; its digest must agree between native, WASM SIMD and forced scalar
WASM (`OCCT_GEOMALGOEXT_SCALAR`). This is finite test evidence, not a proof for
all NaN payloads or non-default native rounding modes. GTests use closed-form
Bernstein derivatives, tensor polynomials, SPD residuals, plane distances and
an independently sampled displacement bound; they are registered in the toolkit.

### Native fitting impact against the pinned stock base

Run `adm/scripts/extensions/spline-fit-benchmark.sh FORK_BUILD STOCK_BASE_BUILD`
after native static Release builds of this revision and stock `b8f597c6`.
The stock-only binary excludes every extension header and call. The fixture supplies
four straight edges of S(u,v)=s(u,v,uv), s=0.1,1,10. The new fit uses 17 sampled
values per edge and the unchanged energy coefficients 1e-6 and 1e-8; stock
`GeomPlate_BuildPlateSurface` receives the same four exact boundary curves,
17 discretization points, G0 and s*1e-7 positional tolerance, with an XY initial
plane. These are different discretizations of the same boundary problem.
Construction/setup and the independent 33x33 quality survey are outside the clock.
Stock measures plate construction, without conversion to a B-spline; the extension
measures the pole solve, without conversion to a surface. Neither includes topology.

| Scale | Stock plate CPU median ms | Spline fit CPU median ms | Stock max grid distance | Fit max grid distance |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 0.7815 | 0.1900 | 2.096e-16 | 2.187e-12 |
| 1 | 0.7210 | 0.1885 | 3.126e-15 | 2.187e-11 |
| 10 | 0.6870 | 0.1085 | 4.213e-14 | 2.187e-10 |

Each row has ten successes and zero refusals after one warmup. GCC 12.2, -O3,
AMD Ryzen 9 7950X3D, TBB off, shared machine, process CPU time. The extension
was faster on these inputs and less accurate than stock; both meet the supplied
tolerance on this finite survey. This is not a global error certificate or a
claim about general patches. Minima/maxima and boundary distances are in
[spline-fit-measurements.json](adm/scripts/extensions/spline-fit-measurements.json).
Native spline GTests pass 6/6. The band harness gives `27232afcfea55306` for all
48 finite matrix configurations in native scalar, WASM SIMD and forced scalar
WASM. That digest establishes these finite factor/solve comparisons only.

### Rim and compaction reference work

Run `adm/scripts/extensions/module-benchmark.sh fit-aux FORK_BUILD STOCK_BASE_BUILD`.
The rim fixture fits a plane displaced by 0.125s using 81 points; stock GeomPlate
receives the same points as G0 constraints. The quality oracle is distance to that
plane on a separate 33x33 grid. The compaction fixture is an affine plane with
inserted quarter knots; stock removes those knots directly, whereas the new
operation additionally certifies a global displacement enclosure. These workloads
are related, not identical. Setup and geometric surveys are outside the clock.

| Scale 1 fixture | Stock CPU ms | New CPU ms | Stock/new max grid distance |
| --- | ---: | ---: | --- |
| parallel-plane rim | 0.2040 | 0.2175 | 3.80e-15 / 1.89e-8 |
| affine knot compaction | 0.0025 | 1.2890 | 2.48e-16 / 3.14e-16 |

Ten successes/zero refusals per row after one warmup, same CPU/compiler settings
as above. The certified compactor is over 500 times slower on this small case;
plain knot removal does not provide its interval certificate. Rim accuracy is
worse than stock here. All three scales and extrema are in
[fit-aux-measurements.json](adm/scripts/extensions/fit-aux-measurements.json).

The native rim/compaction implementations are prepared for review but the app
retains its original implementations: replacing its hot cancellation polling by
an opaque OCCT progress check failed the zero-increase CPU gate. A frozen paired
slope experiment measured +0.665 ns/check (one-sided 95% upper bound +0.806).
This is a boundary performance limitation, not a numerical disagreement. Native
clients may use these APIs; no claim is made that the app consumes them yet.

## Applying modules separately

These modules do not depend on the 17 corrections. The accepted utility commits
remain unchanged. On this pinned base, apply the following commits in order:

| Module | Commits | Dependencies |
| --- | --- | --- |
| Spline fit, portable band, prepared rim/compaction | `bb5616bf41`, `58430e7900`, `dd709ae866` | accepted interval/Bernstein/bounds helpers; existing TKGeomAlgo→TKBRep dependency |
| Exact-section fitting | `a28ec1fd41` | spline-fit batch and accepted Hermite helper |
| Patch chart/interior | `3d1022a70e`, `4273531a63` | spline execution header; accepted Bernstein/interval helpers; existing TKMesh |
| Benchmark documentation | `01bba3caa2`, `f757b82e4c` | corresponding module headers for candidate binaries only |

Commit `dd709ae866` can be cherry-picked before the chart commit: it changes only
the spline basis status boundary. A module cherry-pick includes its own GTests and
FILES registration. If a package already exists, merge the additive FILES entries;
do not replace unrelated toolkit registrations. Native tests above were run on the
combined pinned branch, not on every independent commit combination or upstream tip.
The module benchmark script accepts separate candidate and stock static build trees.
Do not infer independent upstream-master compatibility from the correction branches.

## Reproducing the native measurements

Use Linux, CMake, Ninja, a C++17 compiler, and static Release libraries.
`adm/scripts/extensions/compare.sh BASE CANDIDATE OUTPUT_DIRECTORY` creates
independent worktrees, builds both, and runs the same benchmark source against
each. For this table BASE is `b8f597c6` and CANDIDATE is `14f92b09`.
To measure a correction alone, create a branch at BASE, cherry-pick only that
correction (and its listed dependency), and pass that branch as CANDIDATE.
The exception and null-curve commits are included in the aggregate candidate;
the native benchmarks do not exercise those WebAssembly/null paths.

The benchmark uses `std::clock()` process CPU time, one warm-up and ten measured
runs per case; setup is outside the clock except sphere construction in the mesh
case. Measurements: AMD Ryzen 9 7950X3D, GCC 12.2, `-O3`, TBB off, six CPU container
limit. These are observed medians, not complexity bounds or universal speedups.
The four finite fixtures do not exercise every correction. A build and unchanged
fixture result are not a triggering regression test for an unexercised correction.

Fixtures have independent checks: a 100×80×2 plate with 32 radius-one through holes
has volume `16000−64π` and area `16720+64π`; 20,000 face classification queries
are checked against the disk complement (19,760 inside). Filleting the top hole
edges by radius r=0.2 removes another `32π[2r²+(5/3)r³−(π/2)(1+r)r²]`.
This follows integrating `π[((1+r)−sqrt(r²−t²))²−1]` from t=0 to r
for each of the 32 quarter-circle profiles. The volume check runs outside
the measured interval. Classification checks the aggregate inside count. A radius-ten sphere meshed at
linear deflection 0.01 and angular deflection 0.1 checks mesh vertices against the
sphere. Its 10,176 triangles and filleted volume 15797.13503629154 were identical
between these two builds. These checks do not certify every topological invariant.

| Aggregate CPU median, ms | Stock same base | Corrections | Ratio |
| --- | ---: | ---: | ---: |
| Plate properties | 4.645 | 4.349 | 0.936 |
| Plate classification | 144.9025 | 4.829 | 0.0333 |
| Plate fillet | 11.2375 | 11.1585 | 0.993 |
| Sphere mesh | 34.020 | 26.987 | 0.793 |

Sample minima/maxima and per-correction values are in
[measurements.json](adm/scripts/extensions/measurements.json).
Every standalone same-base correction built and completed all four fixtures.
The analytic fillet-volume check was added after the first timing series, outside
the clock; a later aggregate rerun also passed it. Timed operation bodies are unchanged.
Noise and regressions in individual observations remain visible there. Each entry records
the ten-sample median and extrema, not every individual timing.

| One correction on stock base, CPU median ms | Properties | Classification | Fillet | Mesh |
| --- | ---: | ---: | ---: | ---: |
| Approx_BSplComputeLine-local-endpoint-parameters | 4.528 | 144.438 | 11.102 | 33.727 |
| BOPAlgo_RemoveFeatures-extend-across-feature | 4.468 | 143.753 | 11.858 | 33.534 |
| BRepAlgo_Loop-sort-vertices | 4.522 | 149.160 | 11.270 | 33.138 |
| BRepCheck_Solid-stop-after-enclosed-region | 4.422 | 141.957 | 11.589 | 33.758 |
| BRepFill_Draft-orientation-last-sample | 4.736 | 143.969 | 11.293 | 32.091 |
| BRepGProp_Gauss-used-errors | 4.368 | 147.476 | 11.118 | 32.648 |
| BRepMesh_Delaun-hilbert-insertion | 4.229 | 140.642 | 11.428 | 26.280 |
| BRepMesh_DelaunayDeflectionControlMeshAlgo-node-normals | 4.488 | 140.654 | 11.019 | 33.898 |
| Geom_OffsetSurfaceUtils-cancellable-progress | 4.447 | 144.507 | 11.335 | 33.918 |
| BRepOffset_MakeOffset-cancellable-tolerances | 4.328 | 146.446 | 10.976 | 33.206 |
| BRepOffset_MakeOffset-input-order | 4.363 | 142.315 | 11.386 | 32.689 |
| BRepTopAdaptor_FClass2d-wire-bounds | 4.533 | 5.003 | 11.198 | 32.855 |
| Extrema_ExtPExtS-spline-across | 4.656 | 145.131 | 11.127 | 32.544 |
| IntTools_BeanFaceIntersector-whole-range | 4.953 | 143.904 | 11.118 | 33.328 |
| IntTools_EdgeFace-line-plane-touch | 4.219 | 146.035 | 11.046 | 32.903 |
| IntTools_FaceFace-plane-across-sweep | 4.532 | 151.938 | 10.979 | 32.733 |
| TopOpeBRepBuild_WireEdgeClassifier-box-apart | 4.662 | 144.457 | 10.789 | 35.466 |

## Corrections, derivations and independent application

For each commit: `git switch -c review-N b8f597c6; git cherry-pick COMMIT`.
All are independent except offset tolerance cancellation, which needs the
`Geom_OffsetSurfaceUtils` progress guard first. Existing OCCT tolerances and
limits are retained unless the entry explicitly identifies an added constant.
These are elementary derivations unless a paper reference is named; they are
not claims of formal machine-checked proofs.

Upstream assessment was against official `master` `53706f78a8640c4407bdcc6752df13a270b159f2`.
Branches `upstream/ID` contain one correction each; the loop sorting branch only
adapts include context. Thirteen branches compiled independently. Offset tolerance
cancellation requires the progress guard and is not independently compilable;
`upstream/BRepOffset-cancellable-with-progress` contains both commits and also builds.
Three changes conflict with redesigned upstream implementations and are left
unported rather than making an unsupported semantic adaptation.

### Approx_BSplComputeLine-local-endpoint-parameters

Commit `806f921d8abf23864e4ef0de270f4ac340e653db`. Upstream: applies; native toolkit build passed.

Problem: Endpoint fitting previously mixed physical parameters with normalized local coordinates.

Derivation and limits: For a positive local span h=b-a, s=(t-a)/h and dC/dt=(dC/ds)/h. Solve the endpoint Bezier fit on [0,1] and rescale both endpoint tangents by h. Assumes strictly increasing endpoint parameters. No new threshold or sampling count.

### BOPAlgo_RemoveFeatures-extend-across-feature

Commit `af44c0a7023a37f38c904010e99b5c64630523bb`. Upstream: applies; native toolkit build passed.

Problem: Extending both directions unnecessarily enlarges adjacent face domains during feature removal.

Derivation and limits: On a forward face an oriented parameter tangent (du,dv) has outward transverse direction (dv,-du). Extend only coordinates with a nonzero outward component; a missing pcurve retains the full extension. This local derivation does not prove that five samples capture every tangent on a curved edge. Four subintervals are an unproven sampling choice. Existing angular and geometric resolutions are retained.

### BRepAlgo_Loop-sort-vertices

Commit `2c9b2590d8d4f695a5a33899ed545f334c233978`. Upstream: applies; native toolkit build passed.

Problem: Repeated parameter evaluation and adjacent exchange sorting take quadratic work.

Derivation and limits: Deduplicate in encounter order, evaluate each remaining parameter once, then stable-sort it. For finite parameters this preserves the stable order of the adjacent exchange sort, including equal parameters, while reducing sorting comparisons from quadratic to O(n log n). Assumes the shape map equality agrees with the old duplicate comparison. No new geometric constants.

### BRepCheck_Solid-stop-after-enclosed-region

Commit `561fa613b4d90aa2d0ce3db2c38a53b188b0cb2b`. Upstream: applies; native toolkit build passed.

Problem: Continuing after an enclosed-region error reaches classifiers that were never initialized.

Derivation and limits: Once a second growth shell adds EnclosedRegion, later classifiers have not been initialized. Return before their use. The already-recorded invalidity is retained; no claim is made about collecting further diagnostic errors. No new constants.

### BRepFill_Draft-orientation-last-sample

Commit `2401e8b9212b9695045f34ade8f006c8daa0425e`. Upstream: applies; native toolkit build passed.

Problem: Rounded final abscissae can exceed the wire length and query an absent location law.

Derivation and limits: Exact sample positions satisfy 0 <= (i-1)*length/N <= length. Floating-point multiplication can place the endpoint outside this closed interval. min(sample,length) restores the endpoint bound and prevents querying an absent location law. Interior arithmetic and sample count are unchanged. No new constants.

### BRepGProp_Gauss-used-errors

Commit `bcf74fe6c586fd395899c34da690a9498cdba0ce`. Upstream: not ported: changed upstream implementation.

Problem: Adaptive integration searched unwritten zero-initialized error capacity.

Derivation and limits: The error vector has an explicitly zero-initialized suffix. Its maximum is the prefix maximum unless that is negative, in which case the first suffix zero wins. Scan the used prefix with the same strict comparison and reproduce that suffix case. This retains the selected integration interval and arithmetic. No new numerical constants.

### BRepMesh_Delaun-hilbert-insertion

Commit `3fa3feb2d4fead0288224bc6812f68087c01f4a2`. Upstream: applies; native toolkit build passed.

Problem: Geometrically scattered insertion order makes mesh refinement expensive.

Derivation and limits: Map the parameter box to a 16-bit square grid and order nodes by recursive quadrant rotation, breaking ties by node index. This is an ordering change only; deflection controls are unchanged. Locality is a performance heuristic, not a runtime theorem. The 16-bit grid (65535 maximum, loop starts at bit 15) is an unproven resolution choice; triangle connectivity and sampling may change.

### BRepMesh_DelaunayDeflectionControlMeshAlgo-node-normals

Commit `e80dd7d596030b3dfead2994718d167fea2ba65d`. Upstream: not ported: changed upstream implementation.

Problem: Repeated normal estimation at unchanged mesh nodes duplicates work.

Derivation and limits: For an unchanged surface and identical parameter pair, NormEstim is the same deterministic function call. Cache its value and defined flag keyed by node index plus exact parameter coordinates, invalidating between optimizations. Coordinate validation prevents stale results after index reuse. No tolerance, angular limit or normal formula changes.

### Geom_OffsetSurfaceUtils-cancellable-progress

Commit `6593918e73854db64fd1d3e8729956320d05bca2`. Upstream: applies; native toolkit build passed.

Problem: Offset normal retries can stop making representable progress without observing cancellation.

Derivation and limits: An iteration with unchanged stored parameters evaluates the same normal again and cannot advance. Reject nonfinite or nonrepresentable steps without an iteration bound; retain arithmetic for progressing steps. A borrowed thread-local Message_ProgressRange is scoped by RAII, so nested evaluations restore their callers even on exceptions. No new numerical constants.

### BRepOffset_MakeOffset-cancellable-tolerances

Commit `bdd6b9a822298b88f5792a3ca51bb44c272ab112`. Upstream: applies textually; requires progress-guard dependency to build.

Problem: Edge tolerance sampling lacked cancellation and loaded derivative-dependent adaptors unnecessarily.

Derivation and limits: Positional D0 evaluation uses only the loaded pcurve and surface, not the endpoint derivative adaptors. Load curve then surface and retain all positional samples and arithmetic. Check the progress range before applying incomplete tolerance results. Depends on the offset evaluation progress guard. The new progress weight 1 is an unproven presentation share, not a geometry bound or deadline.

### BRepOffset_MakeOffset-input-order

Commit `b0456807767b5c4e2017882573a2c92cdce0ecfe`. Upstream: applies; native toolkit build passed.

Problem: Offset root traversal depended on allocation-sensitive map order.

Derivation and limits: Sort map entries by TopExp input index before building their generated faces. The traversal becomes a function of input topology instead of allocation addresses, provided every key has an input index. Stable sorting retains existing order for equal indices. No tolerance or geometric constant changes.

### BRepTopAdaptor_FClass2d-wire-bounds

Commit `c8c90a2c33699ec7e532d9614a42ba953fc36357`. Upstream: not ported: changed upstream implementation.

Problem: Using a face-wide box for each wire prevents rejecting unrelated holes cheaply.

Derivation and limits: A polygon lies in the coordinate-wise extrema of its vertices. A query outside that box beyond the classifier tolerance cannot be inside the polygon. Use that wire box for quick rejection while retaining the identical polygon ray test inside and the old fallback for zero-area bounds. No new constants.

### Extrema_ExtPExtS-spline-across

Commit `a71c466051d1aba6fe3aa7227f5bea26488e921e`. Upstream: applies; native toolkit build passed.

Problem: A generic surface extremum solve ignores separability of transverse spline extrusions.

Derivation and limits: For unit D and C(u).D=h constant, S(u,v)=C(u)+vD gives squared distance ||P-h_proj-C(u)||^2+(v-v0)^2. The interior minimum therefore separates into a curve projection and v=v0. Keep the generic path for boundary or unresolved cases. Existing capacity of four extrema is retained: this does not promise every boundary or non-nearest extremum. Existing tolerances are unchanged.

### IntTools_BeanFaceIntersector-whole-range

Commit `dc671c167bc420983b58647b6871781decc0e297`. Upstream: applies; native toolkit build passed.

Problem: Localized edge/face searches may omit intersections outside their selected neighborhood.

Derivation and limits: Use the existing whole-range search for all surface kinds instead of choosing a localized search for certain polynomial surfaces. Both seek the same zero distance locus; no theorem guarantees completeness of the numerical global search. This preserves the previously selected implementation, not a new intersection algorithm. No new constants.

### IntTools_EdgeFace-line-plane-touch

Commit `def5b1740c39f3d443ba016837ee18cf13f35393`. Upstream: applies; native toolkit build passed.

Problem: Numerical touch detection can miss the analytic interior intersection of a line and a plane.

Derivation and limits: For L(t)=p+td and plane (x-q).n=0, t=(q-p).n/(d.n) when d.n is nonzero. Use that intersection only inside both trimmed parameter ranges, retaining the generic extrema search elsewhere. Angular nonparallel recognition uses the existing Precision value. No new thresholds.

### IntTools_FaceFace-plane-across-sweep

Commit `5ce72fca97a8b6313e019f52cb876d884ca96e0f`. Upstream: applies; native toolkit build passed.

Problem: The generic section search loses the exact translated-curve construction of transverse extrusions.

Derivation and limits: For S(u,v)=C(u)+vD, D parallel to plane normal n, and C(u).D constant, v=(q-C(u0)).n/(D.n) is constant. Translate the basis curve by vD; its extrusion pcurve is (u,v). Preserve boundedness for section endpoints, and retain the general path near trimmed ends. Assumes positive-weight spline poles for the constant-height hull argument. Existing OCCT tolerances only.

### TopOpeBRepBuild_WireEdgeClassifier-box-apart

Commit `14f92b09b402c6fcfcbc95129e6a329b4184372e`. Upstream: applies; native toolkit build passed.

Problem: Wire classification performs expensive curve tests even for disjoint parameter boxes.

Derivation and limits: A loop is contained in the bounds of all its existing pcurves. If two expanded boxes are disjoint, one loop is in the unbounded complement of the other box and can reuse its outside classification. Missing pcurves and block classification keep the old path. The outside probe uses xmax+width+1; 1 is an unproven positive displacement, not an error bound. UV expansion retains stored edge tolerances and therefore inherits their unit assumptions.

The Hilbert ordering uses the locality concept discussed in [4]; its grid resolution
remains a heuristic.

## Finite mathematical utilities

The originally extracted utility subset consists of header-only helpers using
OCCT and the C++ standard library. These helpers neither report application
errors nor own an application journal. Their bounded local calculations have no
long-running solver loop needing a progress range. The larger solver modules
described above have separate native failure and progress contracts; their
application integration remains subject to the stated boundary gates.

* `GeomAlgoExt_Bernstein.hxx`: for Bernstein basis `B_i^n`,
  `B_i^n B_j^m = binom(n,i) binom(m,j) / binom(n+m,i+j) B_(i+j)^(n+m)`.
  Degree elevation is the corresponding multiplication by the constant one.
  Tensor-product extension applies the identity in each coordinate [1].
  The cached Pascal table size 96 is an **unproven performance cutoff**.
  Floating coefficient accumulation is not a certified global error bound.
* `GeomAlgoExt_Interval.hxx`: for finite stored IEEE binary64 endpoints, evaluate
  endpoint operations and expand each rounded result by one adjacent representable
  number. The enclosure follows correct round-to-nearest and monotonicity; a
  product needs the four endpoint combinations [2]. Reject nonfinite results.
  This assumes no unsafe reassociation, overflow acceptance, or altered rounding.
  Binary64 precision 53 is derived from the representation, not a fitted constant.
* `GeomAlgoExt_UniformKnots.hxx`: a degree-p clamped spline has endpoint
  multiplicity p+1; its interior simple knots divide the finite interval into
  equal spans. Assumes nonnegative degree, positive span count, finite endpoints.
* `GeomAlgoExt_Hermite.hxx`: convert value and derivative constraints into rows of
  the B-spline collocation matrix using OCCT's basis evaluation, then solve the
  three coordinate right-hand sides. Knot multiplicities and derivative scaling
  follow spline differentiation [3]. A singular matrix returns failure; neither
  invertibility nor floating conditioning is promised for arbitrary constraints.
* `GeomAlgoExt_CubicThrough.hxx`: for distinct increasing sites, a natural C²
  cubic interpolant minimizes integral |s''|² over H² interpolants. Writing
  another interpolant as s+e, integrate each span twice: the fourth derivative is zero, e vanishes at
  sites, s'' is continuous internally and zero at the ends. Thus integral s''e''=0
  and E(s+e)=E(s)+E(e). The implemented tridiagonal system and cubic Bernstein
  controls follow the first/second derivative equations; coefficients 2, 3, 1/2
  and 3/2 and degree/multiplicity 3/4 are derived. Closed and tangent-constrained
  variants have their own endpoint equations, not the natural-boundary theorem.
  `unify` retains an **unproven 1e-12 knot merge threshold**.
* `GeomAlgoExt_ParameterRectangle.hxx`: a positive-weight rational B-spline lies
  in its control hull [1,3]. Bound each parameter pcurve and recognize a full
  rectangle under the existing valid-face and boundary assumptions. This is not
  a global topology certificate. **64 epsilon** and **1e-12 times span** are
  inherited **unproven recognition allowances**, not mathematical bounds.
* `GeomAlgoExt_SurfaceBounds.hxx`: split a B-spline into Bezier patches and take
  coordinate-wise control hull extrema [1,3]. Positive weights are required for
  the rational hull argument. The inherited conversion allowance **1e-10 is
  unproven**; the hull theorem does not prove the conversion error.

Native time for the new helper APIs versus stock has not been measured: stock has
no identical entry points. The integration comparison isolates relocation without
changing their implementation; the first two batches produced identical WASM bytes.
The native utility oracle runner is `adm/scripts/extensions/utilities.sh BUILD_DIRECTORY`.
It checks independent polynomials, long-double inclusions, endpoint jets and hull containment.

No numerical threshold was changed by relocating these utilities. The source
identifiers changed mechanically to the package namespace. Host integration
proof compares exact serialized B-reps and typed refusals against the installed
starting release, rather than rebaselining changed geometry.

## Sources actually consulted

[1] Rida T. Farouki (2012), “The Bernstein polynomial basis: A centennial
retrospective”, Computer Aided Geometric Design 29(6), 379–419.
[Author manuscript](https://www.researchgate.net/profile/James-Peters-3/post/Since_the_Bernstein_basis_is_more_stable_than_the_1_x_x2_x3_basis_why_do_we_use_the_second_one/attachment/59d62e4cc49f478072e9ef35/AS:273573887840257@1442236505726/download/Farouki-on-Berstein-basis-2012.pdf),
degree elevation/product formulas and convex hull sections. The floating-point
implementation does not inherit an exact-arithmetic theorem without rounding analysis.

[2] Siegfried M. Rump (2010), “Verification methods: Rigorous results using
floating-point arithmetic”, Acta Numerica 19, 287–449.
[Author manuscript](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf),
sections 2, 4 and 5.5 on IEEE arithmetic, directed rounding and inclusion.
The finite stored-double helper is narrower than general interval arithmetic.

[3] Carl de Boor (1972), “On calculating with B-splines”, Journal of Approximation
Theory 6(1), 50–62. DOI [10.1016/0021-9045(72)90080-9](https://doi.org/10.1016/0021-9045(72)90080-9).
[Read paper](https://ftp.cs.wisc.edu/Approx/bsplcalc.pdf), recurrence and convex combinations.

[4] John Skilling (2004), “Programming the Hilbert curve”, AIP Conference
Proceedings 707, 381. DOI [10.1063/1.1751381](https://doi.org/10.1063/1.1751381).
[Read transcription](https://doczz.net/doc/7960692/programming-the-hilbert-curve).
This supports the locality construction concept. The implemented two-dimensional
quadrant rotation is not a claim to reproduce Skilling's global Gray-code algorithm;
no speed theorem follows from the citation.

## License and contribution review

New utility files identify Timon Gentzsch as author and use OCCT's LGPL-2.1 terms
with its exception. Existing OCCT files retain their existing notices.
Upstream's own [contribution workflow](dox/contribution/contribution_workflow/contribution_workflow.md)
states: “Before contributing, you must sign the Contributor License Agreement (CLA).”
It requests contributions to IR, a `<Group> - <Summary>` title and review/testing.
No CLA has been signed by preparing this local branch. Ownership and authority to
license the extracted application sources must be confirmed before publication.
Original OCCT project information is preserved below.

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
