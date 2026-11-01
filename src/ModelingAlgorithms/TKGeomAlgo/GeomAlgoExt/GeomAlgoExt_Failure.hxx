// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <TopoDS_Shape.hxx>
#include <stdexcept>
#include <string>
#include <vector>

namespace loft {
// Stable reasons and codes are fields, never tags in diagnostic text.
inline constexpr const char* invalid_input = "invalid_input";             // caller error: wrong shape type, zero vector, ...
inline constexpr const char* unsupported_geometry = "unsupported_geometry"; // valid request the algorithm cannot handle
inline constexpr const char* no_solution = "no_solution";                 // the algorithm ran and produced no usable result
inline constexpr const char* nothing_removed = "nothing_removed";         // defeature left every requested face in place
inline constexpr const char* invalid_result = "invalid_result";           // the result failed the operation's own checks

// Construction diagnostics are separate from the public broad reason and OCCT error code.
// Dispatch must never depend on human-readable text.
enum class FailureKind {
    Construction, NativeBuilder, TangentEdge, UnsupportedRail, InvalidTopology,
    CapConsumed, ClosingConcave, Overflow, SupportContact, Tangency, Curvature,
    RadiusEndpoints, RadiusLaw, MaterialBound
};
const char* failureKindName(FailureKind kind);
int failurePriority(const struct Failure& failure);

/// What `fail` throws. `what()` remains useful for native diagnostics; callers branch on fields.
struct Failure : std::runtime_error {
    std::string operation, reason, code, message;
    std::vector<int> at;
    /// The sub-shapes `at` names, where `fail` was given them: for a caller that names others in
    /// their place (a copy's by the input's) by shape, as the hashes are not unique.
    std::vector<TopoDS_Shape> shapes;
    FailureKind kind;
    Failure(const char* op, const char* why, const std::string& detail, std::vector<int> where,
            const char* stableCode = "CONSTRUCTION_FAILED", FailureKind failureKind = FailureKind::Construction)
        : std::runtime_error(std::string(op) + ": " + detail), operation(op), reason(why),
          code(stableCode), message(detail), at(std::move(where)), kind(failureKind) {}
};

/// The bound of the hashes in Failure::at: the largest, which is the one core asks with.
inline constexpr int failureHashBound = 2147483647;
/// A sub-shape's hash in Failure::at.
int failureHash(const TopoDS_Shape& shape);

[[noreturn]] void fail(const char* operation, const char* reason, const std::string& message);
/// The same, at these sub-shapes of the input.
[[noreturn]] void fail(const char* operation, const char* reason, const std::string& message,
                       const std::vector<TopoDS_Shape>& at, const char* code = "CONSTRUCTION_FAILED",
                       FailureKind kind = FailureKind::Construction);
[[noreturn]] inline void fail(const char* operation, const char* reason, const std::string& message,
                             const std::vector<TopoDS_Shape>& at, FailureKind kind) {
    fail(operation, reason, message, at, "CONSTRUCTION_FAILED", kind);
}

} // namespace loft
