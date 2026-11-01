// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The host seam of the operations in TKAlgoExt. A host keeps shapes under 32-bit handles in a
// store of its own, here called OcctKernel. The algorithms never look inside it: they ask for
// the shape behind a handle, register a result, and release a temporary, through the four
// functions below, which the host defines. The same code therefore runs over any store.

#include <TopoDS_Shape.hxx>

#include <cstdint>
#include <vector>

class OcctKernel;

namespace loft {

/// A result and what became of the tracked input faces: [inputFaceHash, count, resultFaceHash...]
/// records in `modified` and `generated`, and the hashes of input faces in `deleted`.
struct EvolutionData {
    uint32_t resultId = 0;
    std::vector<int> modified;
    std::vector<int> generated;
    std::vector<int> deleted;
};

/// Register a shape in the host's store and return its handle.
uint32_t store(OcctKernel& kernel, const TopoDS_Shape& shape);

/// The shape behind a handle. Throws a typed invalid_input failure with code INVALID_SHAPE_ID.
const TopoDS_Shape& get(const OcctKernel& kernel, uint32_t id);

/// Drop a handle the algorithm itself registered for a temporary.
void release(OcctKernel& kernel, uint32_t id);

/// A result that is a new shape and no face history: `shape`, registered in the store.
inline EvolutionData stored(OcctKernel& kernel, const TopoDS_Shape& shape) {
    EvolutionData out;
    out.resultId = store(kernel, shape);
    return out;
}

} // namespace loft
