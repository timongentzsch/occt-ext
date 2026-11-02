<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# TKAlgoExt

Modelling operations built on OCCT's topology, geometry and Boolean toolkits. Nine packages, all in namespace
`BRepExt`:

| Package | Operations |
| --- | --- |
| `BRepBlendExt` | Rolling-ball blend solver and its attachment to the body |
| `BRepFilletExt` | Constant, variable, range and station fillets, chamfer, blend validity checks |
| `BRepPatchExt` | Boundary-first fill of a loop of edges with G0 / G1 / G2 to neighbouring faces |
| `BRepSurfaceExt` | Alignment, bridging, deformation, face separation proofs |
| `BRepCurveExt` | Offset, outline, corner rounding, G0 to G3 bridge curves, projection, regions |
| `BRepOffsetExt` | Solid offset, hollow, shell, thicken, push-pull, defeaturing |
| `BRepSkinExt` | Gordon and guided skins, framed sweeps, draft |
| `BRepBooleanExt` | Cell-based Booleans with checked fallbacks |
| `BRepToolsExt` | Properties, boxes, dimensions, hidden line, shared checks and history |

Construction, mathematics, sources and unproven choices of each are in
[dox/extensions/operations.md](../../../dox/extensions/operations.md). The interface to an embedding program
(`BRepToolsExt_Session.hxx`) and the build options are in the top-level [README](../../../README.md).
Tests are in `GTests/`, one file per package.
