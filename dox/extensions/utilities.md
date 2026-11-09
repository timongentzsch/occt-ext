<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# Numerical utilities, fitting and Patch charts

Two toolkits gain packages: `GeomAlgoExt` in `TKGeomAlgo` (geometry and numerics) and `BRepAlgoExt` in
`TKMesh` (charts, trim integration, properties). Both use existing OCCT toolkits only and are independent of
`TKAlgoExt`, except that the blend and Patch solvers call them. Each section states the mathematics, the
assumptions under which it holds, and the numbers that are choices rather than theorems. "Unproven" always
means a number or rule chosen for the code, not derived. Sources are listed in [references.md](references.md).

## Boundary-first patch chart

`BRepAlgoExt_PatchFlow` builds a chart of a disc-like boundary loop. It uses `TKMesh` because the chart needs its
Delaunay mesher. The interior survey used by Patch is `BRepPatchExt_Inside`; its ratio-four acceptance rule is an
unproven heuristic ([BRepPatchExt.md](operations/BRepPatchExt.md)).

The cotangent stiffness matrix L is formed, the interior equation L_II u = K_I is solved with zero boundary
log scale, and the boundary turns follow from pi - angleSum - (Lu)_boundary. For boundary directions t_i and
original lengths l_i, the closure minimises sum (delta l_i)^2 / l_i subject to sum (l_i + delta l_i) t_i = 0.
Lagrange multipliers give delta l_i = -l_i t_i^T G^-1 sum l_j t_j with G = sum l_i t_i t_i^T, which is the
implemented 2x2 closure. A harmonic x/y extension then solves the Dirichlet energy equations. This follows
Sawhney and Crane (2017), *Boundary First Flattening*, sections 3 to 5, especially equations 6, 14 and 20. Their
section 4.4 separates boundary interpolation from harmonic conjugation; this implementation interpolates, so
**exact discrete conformality is not guaranteed**.

Assumptions: a finite regular surface, a connected non-degenerate disc mesh, a usable boundary and non-singular
closure, positive resulting edge lengths and triangle orientation. Cotangent Dirichlet systems are positive
definite with a fixed non-empty boundary under the usual finite-element assumptions; arbitrary inputs can
refuse. The node cap, curvature survey and parameter-space lattice do **not** certify 3D deflection or global
injectivity. The chord estimate k l^2 / 8 is local and asymptotic, not a global surface error bound. The chart
may extrapolate off its mesh.

Unproven choices: 32 000 chart nodes; conjugate-gradient residual target 1e-13 and 4n passes; at least two edge
subdivisions, a half-spacing rim margin, and sqrt(triangleCount / 2) for the lookup grid. The factors sqrt(3)/2 and
sqrt(3)/4 of the triangular lattice are derived from equilateral geometry; the cotangent 1/2 and triangle 1/2
are derived energy and area factors. The interior survey has its own sample counts and guards (see its source).

A native test checks planar pair distances at three scales. Reproduce the comparison with
`adm/scripts/extensions/module-benchmark.sh patch-chart FORK_BUILD STOCK_BASE_BUILD` after building static
Release `TKMesh` from both revisions. The fixtures are a rectangle and a quarter cylinder at scales 0.1, 1 and
10. Stock offers meshing and scaled carrier UVs but no boundary-first chart, so the reference does less work and
produces different meshes. Geometry and angle checks are outside the CPU clock.

| Scale 1 fixture | Stock CPU ms | New CPU ms | Stock / new max cosine error | Stock / new nodes |
| --- | ---: | ---: | --- | --- |
| plane-chart | 0.0760 | 0.2770 | 0 / 5.54e-13 | 4 / 75 |
| cylinder-chart | 0.1005 | 0.4510 | 0.000155 / 0.00222 | 18 / 101 |

Ten successes and no refusals per row after one warm-up; GCC 12.2 `-O3`, TBB off, AMD Ryzen 9 7950X3D, shared
machine. The chart can be slower than meshing with carrier UVs alone, and this does not compare equivalent
flattening solvers. All scales, extrema, relative edge distortions and curvature errors are in
[patch-chart-measurements.json](../../adm/scripts/extensions/patch-chart-measurements.json).

## Constrained tensor-spline fitting and certified compaction

`GeomAlgoExt_PatchFit`, `GeomAlgoExt_PatchRim` and `GeomAlgoExt_PatchCompact` extend `TKGeomAlgo`;
`GeomAlgoExt_Band` is their banded Cholesky engine. They accept OCCT geometry and ordinary C++ containers.
`BRepExt::Failure`, `fail`, `FailureKind` and the stable reason constants are in `GeomAlgoExt_Failure`.
`Failure` derives from `std::runtime_error`; substituting `Standard_Failure` would change which existing OCCT
`catch` blocks intercept a refusal, so adapting that taxonomy is a task for upstream review. The integer
diagnostic ranking is categorical policy, not geometry; its scores 100/90/80/70/60/40/30/10/0 encode a priority
order. The hash modulus 2147483647 is the signed 31-bit interface ceiling.

Cancellation is a direct call of `sessionCancelled()` at every polling site. A native caller constructs
`BRepExt::ProgressScope(range)` from `GeomAlgoExt_Progress.hxx` around the entry call; its thread-local range
backs the default hook, and nested scopes restore the previous range. `BRepExt_PatchPhase` is a separate phase
observer hook whose native default is a no-op. No polling or erased callback is inserted into a hot loop. The
native default implementation is compiled only for native builds: in a WebAssembly archive a native definition
would satisfy the linker before an embedder's import and silently replace cancellation and observation.

*Fit.* For a tensor spline S(u,v) = sum P_ij N_i(u) M_j(v), scalar conditions are linear forms a_r^T p = b_r in
the free pole coordinates. The energy is

    integral [ alpha (|S_uu|^2 + 2|S_uv|^2 + |S_vv|^2) + beta (|S_uuu|^2 + 3|S_uuv|^2 + 3|S_uvv|^2 + |S_vvv|^2) ] du dv.

The coefficients 1, 2, 1 and 1, 3, 3, 1 count ordered mixed derivatives, so they are derived. Substituting the
basis gives a positive semidefinite quadratic p^T H p + 2 g^T p + c. Weighted constraint rows add A^T W A;
eliminating fixed poles adds their exact contribution to the right-hand side. Cholesky solves the resulting
symmetric positive definite system, and multipliers correct row residuals using products with A M^-1 A^T.
Stopped iterations are returned as they stand; convergence, consistency and global optimality are not promised.
The semidefinite claim assumes non-negative energy coefficients, row weights and density, finite regular knot
spans and consistent constraint dimensions; solvability requires positive definiteness after fixed-pole
elimination. These are preconditions of the caller, not a proof that arbitrary inputs give a usable surface.
Gauss quadrature with p + 1 points integrates the polynomial basis products on each knot span exactly;
arbitrary supplied density functions have no such guarantee. Basis locality and derivative identities follow
de Boor (1972), sections 1 and 2. These derivations do not certify conditioning.

*Rim.* The correction linearises the nearest-point displacement along a supplied unit vector and solves
min_x ||A x - r||^2 + d^2 ||x||^2, whose normal equations are (A^T A + d^2 I) x = A^T r; unobserved poles stay
fixed.

*Compaction.* OCCT knot removals are proposed, the original and the candidate are converted to a common piecewise
Bernstein basis with outward-rounded intervals, and the displacement is bounded by max_i ||Delta P_i||_2.
Non-negative basis functions that sum to one imply this bound over the whole rectangle, provided both conversions
are enclosed, degrees and full knot sequences agree, weights are non-rational and the input is finite, clamped
and non-periodic. No normal, trim or continuity bound follows. The outward interval enclosure uses the inclusion
principle in Rump (2010), section 5.5; the reference does not certify the implementation outside its stated
rounding and finite-input assumptions.

*Choices that are not theorems.* Fit diagonal hold 1e-13, factor 100 and cap 1e-6; stall count 12, stall share
0.02, at most 400 passes, energy-ratio stop 1e-20 (not binary64 epsilon squared); basis capacity of degree 5 and
five terms per condition. Rim: 8 Newton steps, determinant cutoff 1e-12, squared-step cutoff 1e-28, sliver ratio
0.2, three solves, damping 0.03, slow-progress ratio 0.5, cosine cutoff 0.17. Compaction: 512 trials, 16 removals,
262 144 coefficients, 8 x 1024 x 1024 visits, a proposal multiplier of 1024, a 33-point rejection screen and a
1e-10 screen allowance. Cancellation is polled every 1024 rows in rim fitting and every 16 or 256 items in
compaction; these are scheduling choices, not latency bounds. The interior acceptance heuristic and the flat G2
limit of Patch remain outside these modules and unproven.

*Portable SIMD.* The WebAssembly SIMD path is kept. A portable two-lane path rounds each operation through
volatile binary64 stores to prevent FMA contraction. The harness `adm/scripts/extensions/band-parity.cxx`
exercises paired and ordinary factors and solves; its digest must agree between native, WebAssembly SIMD and
forced scalar WebAssembly (`OCCT_GEOMALGOEXT_SCALAR`). The digest `27232afcfea55306` is the same for all 48
finite matrix configurations in all three modes. This is finite test evidence, not a proof for all NaN payloads
or non-default rounding modes. The GTests use closed-form Bernstein derivatives, tensor polynomials, SPD
residuals, plane distances and an independently sampled displacement bound.

### Fitting against stock

Run `adm/scripts/extensions/spline-fit-benchmark.sh FORK_BUILD STOCK_BASE_BUILD` after static Release builds
of this revision and of stock `b8f597c6`. The stock binary excludes every extension header. The fixture supplies
four straight edges of S(u,v) = s (u, v, uv) for s = 0.1, 1, 10. The fit uses 17 sampled values per edge and the
energy coefficients 1e-6 and 1e-8; stock `GeomPlate_BuildPlateSurface` receives the same four boundary curves,
17 discretisation points, G0 and a positional tolerance of s * 1e-7, with an XY initial plane. These are
different discretisations of one boundary problem. Setup and the independent 33 x 33 quality survey are outside
the clock. Stock is timed for plate construction without conversion to a B-spline; the extension for the pole
solve without conversion to a surface. Neither includes topology.

| Scale | Stock plate CPU median ms | Spline fit CPU median ms | Stock max grid distance | Fit max grid distance |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 0.7815 | 0.1900 | 2.096e-16 | 2.187e-12 |
| 1 | 0.7210 | 0.1885 | 3.126e-15 | 2.187e-11 |
| 10 | 0.6870 | 0.1085 | 4.213e-14 | 2.187e-10 |

Each row has ten successes and no refusals after one warm-up (GCC 12.2 `-O3`, Ryzen 9 7950X3D, TBB off, process
CPU time). On these inputs the fit is faster and less accurate than stock; both meet the supplied tolerance.
This is not a global error certificate or a claim about general patches. Extrema and boundary distances are in
[spline-fit-measurements.json](../../adm/scripts/extensions/spline-fit-measurements.json). The native spline
GTests pass 6 of 6.

Rim and compaction are measured by `module-benchmark.sh fit-aux FORK_BUILD STOCK_BASE_BUILD`. The rim fixture
fits a plane displaced by 0.125 s through 81 points; stock `GeomPlate` receives the same points as G0
constraints and the oracle is the distance to that plane on a separate 33 x 33 grid. The compaction fixture is an
affine plane with inserted quarter knots; stock removes those knots directly, whereas the new operation also
certifies a global displacement enclosure. The workloads are related, not identical.

| Scale 1 fixture | Stock CPU ms | New CPU ms | Stock / new max grid distance |
| --- | ---: | ---: | --- |
| parallel-plane rim | 0.2040 | 0.2175 | 3.80e-15 / 1.89e-8 |
| affine knot compaction | 0.0025 | 1.2890 | 2.48e-16 / 3.14e-16 |

The certified compactor is over 500 times slower on this small case; plain knot removal does not give its
interval certificate. Rim accuracy is lower than stock's here. All scales and extrema are in
[fit-aux-measurements.json](../../adm/scripts/extensions/fit-aux-measurements.json).

## Finite mathematical utilities

Header-only helpers using OCCT and the C++ standard library. Their bounded local calculations need no
progress range. The native oracle runner is `adm/scripts/extensions/utilities.sh BUILD_DIRECTORY`; it checks
independent polynomials, long-double inclusions, endpoint jets and hull containment. Stock has no entry points
equivalent to these helpers, so no timing against stock is given.

* `GeomAlgoExt_Bernstein.hxx`: for the Bernstein basis `B_i^n`,
  `B_i^n B_j^m = binom(n,i) binom(m,j) / binom(n+m,i+j) B_(i+j)^(n+m)`. Degree elevation is the
  multiplication by the constant one, and the tensor-product extension applies the identity in each
  coordinate [1]. The cached Pascal table size 96 is an **unproven performance cutoff**. Floating-point
  coefficient accumulation is not a certified global error bound.
* `GeomAlgoExt_Interval.hxx`: for finite stored IEEE binary64 endpoints, each endpoint operation is evaluated
  and its rounded result widened by one adjacent representable number. The enclosure follows from correct
  round-to-nearest and monotonicity; a product needs the four endpoint combinations [2]. Non-finite results are
  rejected. This assumes no unsafe reassociation, overflow acceptance or altered rounding. The precision 53 is
  derived from the representation.
* `GeomAlgoExt_UniformKnots.hxx`: a degree-p clamped spline has end multiplicity p + 1, and its interior simple
  knots divide the finite interval into equal spans. Assumes non-negative degree, positive span count and finite
  endpoints.
* `GeomAlgoExt_Hermite.hxx`: value and derivative constraints become rows of the B-spline collocation matrix,
  using OCCT's basis evaluation, and the three coordinate right-hand sides are solved. Knot multiplicities and
  derivative scaling follow spline differentiation [3]. A singular matrix returns failure; neither invertibility
  nor conditioning is promised for arbitrary constraints.
* `GeomAlgoExt_CubicThrough.hxx`: for distinct increasing sites, the natural C2 cubic interpolant minimises
  integral |s''|^2 over H2 interpolants. Writing another interpolant as s + e and integrating by parts twice (the
  fourth derivative is zero, e vanishes at the sites, s'' is continuous in the interior and zero at the ends)
  gives integral s'' e'' = 0 and so E(s + e) = E(s) + E(e). The tridiagonal system and the cubic Bernstein
  controls follow the first- and second-derivative equations; the coefficients 2, 3, 1/2 and 3/2 and the
  degree and multiplicity 3 and 4 are derived. Closed and tangent-constrained variants have their own end
  equations and not the natural-boundary theorem. `unify` has an **unproven 1e-12 knot merge threshold**.
* `GeomAlgoExt_ParameterRectangle.hxx`: a positive-weight rational B-spline lies in its control hull [1, 3].
  Each parameter pcurve is bounded and a full rectangle recognised under the usual valid-face and boundary
  assumptions; this is not a global topology certificate. **64 epsilon** and **1e-12 times the span** are
  recognition allowances that are not mathematical bounds.
* `GeomAlgoExt_SurfaceBounds.hxx`: a B-spline is split into Bezier patches and the coordinate-wise control hull
  extrema are taken [1, 3]. Positive weights are required for the rational hull argument. The conversion
  allowance **1e-10 is unproven**; the hull theorem does not bound the conversion error.

## Surface differential jets and fit diagnostics

`GeomAlgoExt_SurfaceJet` computes first and second derivatives, the oriented normal, one-sided spline
evaluation and curvature comparison. With J = [Su Sv], I = J^T J and II_ij = n . S_ij, the coefficients of an
ambient tangent t are I^-1 J^T t, so II(a, b) is evaluated in parameter coordinates without assuming an
orthonormal parametrisation. The maximum absolute directional curvature difference of a symmetric 2x2 matrix D is
|tr(D) / 2| + hypot((D11 - D22) / 2, D12). One-sided U derivatives use adjacent distinct knots; V uses its
containing span. Assumes regular differentials and compatible tangent planes. Unproven cutoffs: det(I) <= 1e-24,
normal magnitude < 1e-15 and curvature scale > 1e-15; they are not a scale-invariant singularity theorem. The
polynomial GTest has independent piecewise graph derivatives at C1 joins, including the V-span choice and
reversed normals. The header needs `TKGeomAlgo` and `Failure` only.

`GeomAlgoExt_FitDiagnostics` uses those jets to report fixed-parameter residuals and normal angles, a signed
projected Jacobian, and a composite-trapezoid bending energy (||Suu||^2 + 2||Suv||^2 + ||Svv||^2 in normalised
parameters, scaled by length^2). The 32 x 32 cell grid (1089 vertices), the reference-area squared cutoff 1e-24
and the minimum initialisation 1e100 are operational choices without proof. Sampling is not a global regularity or
residual certificate.

## Planar curve containment

`GeomAlgoExt_PlanarCurve` transforms the trimmed curve into the supplied plane frame and asks OCCT's optimal
curve bound for [zmin, zmax]. The test max(|zmin|, |zmax|) <= tolerance is necessary and sufficient for this
bound to lie in the tolerance slab. It uses the caller's tolerance and OCCT's extrema algorithm, with no new
sampling count or constant, and assumes OCCT's bound encloses the curve; unsupported extrema keep stock
behaviour. It depends on `TKGeomAlgo` and `Failure` only.

## Contour seeds and prism diagnostics

`BRepAlgoExt_Contours` enumerates a maker's propagated contour and chooses a periodic spline's closing knot as a
deterministic seed; it does not change the curve. `BRepAlgoExt_Prism` supplies conservative diagnostics after a
failed construction. A planar set contained in a strip of width w cannot contain a disk of diameter 2r > w; this
proves the cap-collapse implication for the recognised right prisms. The curvature diagnostic solves
det(II - k I) = 0 from exact surface derivatives and is a sampled diagnostic, not a certificate between
samples. Assumes regular support derivatives and correctly recognised prism supports. Unproven allowances:
fineAngle 1e-12; strip slack 2e-7; squared differential cutoffs 1e-20 and 1e-24; classifier steps 1e-4, 1e-3,
1e-2 and tolerance 1e-8; 4096 samples (one for constant-curvature elementary supports); curvature-ratio allowance
1e-7. Both depend on `TKMesh` and `TKGeomAlgo`.

## Trim integration, read-back and stored tolerances

`BRepAlgoExt_TrimDomain` partitions rational boundary arcs into v-monotone pieces. For Y/W, v' has numerator
Y'W - YW'; Bernstein hull subdivision isolates its zeros. Horizontal crossings define the trimmed intervals,
which are integrated with embedded Gauss 7 / Kronrod 15 rules for area and first moments. Volume uses the
divergence identity V = (1/3) integral (p - origin) . n dA. Near square-root end behaviour, the substitution
v = lo + width sin^2(t) regularises the endpoint. Assumes a readable valid trim, a regular finite carrier and
positive weights for the convex-hull argument; classification may refuse as Unknown. The embedded difference is
an error estimate, **not a rigorous enclosure**. The nodes and weights are those of QUADPACK `dqk15`
(Piessens and de Doncker, 1980, revised 1983; table by L. W. Fullerton, 1981), read at
<https://netlib.org/quadpack/dqk15.f>; its affine rule is used, its separate resasc and machine-error
estimator is not. Bernstein identities follow de Boor [3].

Unproven choices: station merge and horizontal-piece cutoff 1e-13; relative root floor 1e-14; root depth 42;
boundary-box depth 40; classifier 1e-12; crossing residual 1e-15, parameter width 1e-16 and 100 iterations;
quadrature depth 18, a two-million call budget, relative budgets 1e-8 (area) and 1e-11 (volume), a 0.1 / 0.9
inner / outer split; scale floor 1e-300.

`BRepAlgoExt_ReadBack` normalises geometry and removes trimming. It searches for fixed points of the stock
reader's normalisation, preserving unchanged input geometry in `settleNew`. Its 8 normalisation iterations,
8-ULP coordinate search, 256 frame rotations at 1e-13 radians and 8 conic settling iterations are search limits
without proof. Existence and global error bounds are not claimed; the fallback is unchanged.

`BRepAlgoExt_Tolerance` returns stored vertex, edge and face tolerances and their maximum over a shape;
containers have zero and non-finite values map to infinity. It introduces no tolerance floor. Their native tests
use exact eccentric-hole moments, nested stored maxima and independent repeated reader normalisation.

`BRepAlgoExt_Properties` is the adaptive boundary property core and the analytic area of rectangular faces of
elementary surfaces. Areas integrate ||Su x Sv||; signed volume integrates (p - origin) . (Su x Sv) / 3.
Rectangles of elementary surfaces use their analytic Jacobian integrals, excluding cone-apex crossings and
spindle tori; other faces use the quadrature and fallback policy. The quadrature is a numerical estimate, not a
rigorous certificate. Unproven choices: Gauss order 8; adaptive interval cutoff 1e-14; noise floors 1e-12 and
1e-14; area relative budget 1e-8; volume relative budget 1e-11 and floor 1e-24; relative volume noise 1e-13 and
noise share 1e-3; inner / outer divisors 10 and 100; derivative cutoff 1e-18; stock area fallback tolerance
1e-9; adaptive depth 18 and two million evaluations. Gauss nodes and weights derive from Legendre polynomial
exactness. The sphere-rectangle oracle uses r^2 du (sin v1 - sin v0).
