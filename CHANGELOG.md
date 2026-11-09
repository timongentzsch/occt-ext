<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# Changelog

All notable changes to the extensions on top of OCCT. The base is OCCT `b8f597c6` (30 July 2026).
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.1.1] - 2026-11-09

### Added

- `BRepBlendExt`: a round whose centre path bends tighter than its radius is trimmed where it crosses itself
  between its contacts.
- A GitHub Actions workflow that builds the extension toolkits and runs their tests on Linux and macOS
  (`.github/workflows/extensions.yml`), with a status badge.
- A Quick start in the README, a changelog, `dox/extensions/limits.md` for the open work, and a README for
  `adm/scripts/extensions`.

### Changed

- The operations reference is one page per package (`dox/extensions/operations/`), with an index.
- Source files of the extension packages carry one header: a copyright line and the SPDX identifier
  `LGPL-2.1-only WITH OCCT-exception-1.0`. The offer of Open CASCADE's commercial licence is dropped from files
  authored here; OCCT's own files keep their notices.
- The upstream all-platform workflow runs by hand only.
- The benchmark sources and shell drivers are formatted and commented.
- Comments and documentation no longer compare with unnamed external modellers.
- The public operations (`BRepToolsExt_Operations.hxx`, `BRepSurfaceExt_Tools.hxx` and their neighbours) name
  every parameter and document each function in a line. Integer mode codes are scoped enumerations with the same
  values (`OffsetJoin`, `CurveExtension`, `PatchMethod`, `PatchFlow` and others); Booleans take `Select` and blends
  `blending::Law`. Source change for callers; no behaviour change.
- `Rung` is `BooleanAttempt` (`BooleanResult::attempts`), and `BooleanRequest::journal` is `subShapeHistory`.
- `BRepBlendExt`: the attachment is split into one source file per stage with a shared header, the solver into
  its entry and its stages, and its longest functions into named steps. No behaviour change.
- `BRepPatchExt`, `BRepSkinExt`, `BRepSurfaceExt`, `GeomAlgoExt`: the Patch fill is split into one source file per
  stage, no function is longer than 150 lines, and the fill's inline tolerances are named constants. No behaviour
  change.

### Fixed

- `BRepOffsetExt`: a kept offset's edge tolerance is the largest measured distance between the edge's curves where
  that stays within 1e-6, so valid offsets are no longer rejected as invalid shapes.
- The package README of `TKAlgoExt` no longer lists the removed bridge curves and corner rounding.

## [0.1.0] - 2026-11-09

First release.

- `TKAlgoExt`, nine packages: rolling-ball fillets and chamfers (`BRepBlendExt`, `BRepFilletExt`), boundary-first
  Patch fills (`BRepPatchExt`), curve and surface operations (`BRepCurveExt`, `BRepSurfaceExt`), offsets, shells and
  thickening (`BRepOffsetExt`), Gordon and guided skins, framed sweeps and drafts (`BRepSkinExt`), cell-based Booleans
  (`BRepBooleanExt`), and shared properties and checks (`BRepToolsExt`).
- Numerical utilities in `TKGeomAlgo` (`GeomAlgoExt`) and `TKMesh` (`BRepAlgoExt`): interval enclosures, Bernstein
  products, constrained tensor-spline fitting, certified knot compaction, trim integration, surface jets.
- Corrections to existing OCCT algorithms, one commit each, with the derivation and limits of each
  ([corrections.md](dox/extensions/corrections.md)).
- Benchmarks against stock OCCT ([benchmarks.md](dox/extensions/benchmarks.md)) and 153 tests with closed-form oracles.
