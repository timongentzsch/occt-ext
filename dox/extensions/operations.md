<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# Operations in TKAlgoExt

One page per package. Each entry gives the construction, the mathematics it rests on, the sources read for it, what is
not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule
chosen for the code, not derived; the constants are named in the source next to their use.
The whole-operation CPU and quality comparison with stock OCCT is in [benchmarks.md](benchmarks.md); the
utilities and their error budgets are in [utilities.md](utilities.md), the literature in
[references.md](references.md), and current limits in [limits.md](limits.md).

| Package | Operations | Page |
| --- | --- | --- |
| `BRepBlendExt`, `BRepFilletExt` | Rolling-ball blend solver, constant, variable, range and station fillets, chamfer, corners, blend validity checks | [BRepBlendExt.md](operations/BRepBlendExt.md) |
| `BRepPatchExt` | Boundary-first fill of a loop of edges with G0 / G1 / G2 to neighbouring faces | [BRepPatchExt.md](operations/BRepPatchExt.md) |
| `BRepCurveExt` | Planar and 3D curve offset, outline, tangent circles, projection, rebuild, refinement, regions, spiral | [BRepCurveExt.md](operations/BRepCurveExt.md) |
| `BRepSurfaceExt` | Alignment, bridging, constraint jets, deformation, face support and rebuild, face separation proofs | [BRepSurfaceExt.md](operations/BRepSurfaceExt.md) |
| `BRepOffsetExt` | Solid offset, hollow, shell, thicken, push-pull, face offset that keeps every face, defeaturing | [BRepOffsetExt.md](operations/BRepOffsetExt.md) |
| `BRepSkinExt` | Gordon and guided skins, framed sweeps, draft | [BRepSkinExt.md](operations/BRepSkinExt.md) |
| `BRepBooleanExt` | Cell-based Booleans with a ladder of fallbacks, checked solids, imprints | [BRepBooleanExt.md](operations/BRepBooleanExt.md) |
| `BRepToolsExt` | Properties, boxes, dimensions, hidden line, shared checks and history; the host interface | [BRepToolsExt.md](operations/BRepToolsExt.md) |
