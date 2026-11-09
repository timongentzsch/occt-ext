<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# Benchmark and verification scripts

Every benchmark is a single C++ file that is compiled against the static libraries of a native Release build
(`libTK*.a` found in the build directory) by its shell driver. CPU time is `std::clock()`: one warm-up and ten
measured runs, reported as median, minimum and maximum. The method and the caveats are in
[benchmarks.md](../../../dox/extensions/benchmarks.md) and [utilities.md](../../../dox/extensions/utilities.md).

## Scripts

| File | What it does |
| --- | --- |
| `operation-benchmark.sh`, `operation-benchmark.cxx` | Whole operations of `TKAlgoExt` against stock OCCT and stock followed by the same validity check, on inputs with a closed-form answer. One row per operation, input and implementation. |
| `operation-chart.py` | Draws the output of `operation-benchmark.sh` as the grouped bar chart `dox/extensions/operation-benchmark.svg` (SVG, no dependencies). |
| `benchmark.sh`, `benchmark.cxx` | CPU time of stock operations that the corrections speed up or fix (properties, face classification, fillet, mesh) on the benchmark plate and sphere. Prints one JSON line per case. |
| `compare.sh` | Builds a base and a candidate revision in two git worktrees and runs `benchmark.sh` on each, to measure a correction against the commit before it. |
| `module-benchmark.sh` | Runs `FAMILY-benchmark.cxx` (family `spline-fit`, `patch-chart` or `fit-aux`) against this tree and against a stock build. |
| `spline-fit-benchmark.sh`, `spline-fit-benchmark.cxx` | Constrained tensor-spline fit (`GeomAlgoExt_PatchFit`) against `GeomPlate`, on boundary nets of known surfaces. |
| `patch-chart-benchmark.cxx` | Boundary-first chart of a loop (`BRepAlgoExt_PatchFlow`) against `BRepMesh` with carrier parameters. Run through `module-benchmark.sh patch-chart`. |
| `fit-aux-benchmark.cxx` | Rim fit and certified knot compaction (`GeomAlgoExt_PatchRim`, `GeomAlgoExt_PatchCompact`). Run through `module-benchmark.sh fit-aux`. |
| `utilities.sh`, `utilities.cxx` | Oracle checks of the numerical utilities (Bernstein, interval, Hermite, cubic, surface bounds) against independent polynomials and long-double inclusions. No timing. |
| `band-parity.cxx` | Digest of the banded solver's factors and solutions. It must agree between native, WebAssembly SIMD and forced scalar builds. |

The benchmark sources compile against a stock tree too, with `-DOCCT_EXTENSION_STOCK_ONLY`, which leaves out
the calls into this tree.

## Recorded results

| File | Content |
| --- | --- |
| `operation-benchmark-results.txt` | Output of `operation-benchmark.sh`: median [min, max] CPU time and the quality measure per row. The source of the chart and of the table in the top-level README. |
| `measurements.json` | Output of `benchmark.sh` and the per-correction runs of the corrected algorithms: CPU time per case and the result value that must agree with stock. |
| `spline-fit-measurements.json` | Spline fit against `GeomPlate`: CPU time, maximum grid error and maximum boundary error per case, scale and method. |
| `patch-chart-measurements.json` | Boundary-first chart: CPU time, triangle count, maximum cosine difference, edge-length error and curvature error per case. |
| `fit-aux-measurements.json` | Rim fit and knot compaction: CPU time and the error measures of each case, with the base commit and clock. |
