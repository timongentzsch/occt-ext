<!-- Copyright (c) 2026 Timon Gentzsch. LGPL-2.1-only WITH OCCT-exception-1.0. -->

# BRepToolsExt: shared host interface

Part of the [operations reference](../operations.md). Each entry gives the construction, the mathematics it rests on, the sources read, what is not proven, and the stock OCCT facility it relates to. "Unproven" always means a number or rule chosen for the code, not derived; the constants are named in the source next to their use.

## Contents

- [Native test host reference lifetime](#native-test-host-reference-lifetime)

## Notes by topic

### Native test host reference lifetime

The native test host keeps shape entries in `std::deque`. Its `get()` interface returns a
reference; an operation can retain that reference while `store()` appends an intermediate
shape. A vector can invalidate it when it reallocates, making an otherwise identical operation
depend on the number of earlier stores. This is a host-lifetime correction, not a change to
the geometric construction. ISO C++ working draft, *Containers library*,
[deque.modifiers], paragraph 1 (accessed 2026), guarantees that end insertion preserves
references to existing elements: <https://eel.is/c++draft/deque.modifiers#1>.
