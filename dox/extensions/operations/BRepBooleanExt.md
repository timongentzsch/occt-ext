<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepBooleanExt: cell-based Booleans

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Booleans](#booleans)
- [Unproven choices](#unproven-choices)

## Notes by topic

### Booleans

Cell selection from one General Fuse with a ladder of fallbacks, each rung accepted only if checked; a sheet tool
divides rather than subtracts.

### Unproven choices

Seam and pairing conventions of closed skins (a convention, not derived), tolerance ladders of the Boolean
fallbacks, the sampling counts of the offset and skin validity checks.
