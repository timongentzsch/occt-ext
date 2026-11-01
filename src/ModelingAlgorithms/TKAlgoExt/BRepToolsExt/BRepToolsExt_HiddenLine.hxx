// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepToolsExt_Operations.hxx>

namespace loft {
// A projected curve retains its source body and line category. Coordinates are in the view's
// plane; only their approximation to a polyline depends on the drawing tolerance.
struct DrawingCurve {
    int body, kind;
    std::vector<double> points;
};
std::vector<DrawingCurve> hiddenLine(OcctKernel &, const std::vector<uint32_t> &, Vec3 normal,
                                    Vec3 x, const std::vector<double> &options);
}
