// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_History.hxx>

#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepTools.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <gp_GTrsf.hxx>

#include <cmath>

namespace BRepExt {

EvolutionData scaleAxes(OcctKernel& kernel, uint32_t shapeId, Vec3 factors, Vec3 center,
                        const std::vector<int>& inputFaceHashes, int hashUpperBound) {
    const char* const op = "scaleAxes";
    const double scale[] = {factors.x, factors.y, factors.z};
    const double origin[] = {center.x, center.y, center.z};
    gp_GTrsf transform;
    for (int i = 0; i < 3; i++) {
        const double translation = (1.0 - scale[i]) * origin[i];
        if (!std::isfinite(scale[i]) || scale[i] <= 0.0 || !std::isfinite(origin[i]) ||
            !std::isfinite(translation)) {
            fail(op, invalid_input, "scale factors must be finite and positive, and the center finite");
        }
        transform.SetValue(i + 1, i + 1, scale[i]);
        transform.SetValue(i + 1, 4, translation);
    }
    FaceHistory history(op, inputFaceHashes, hashUpperBound);
    const TopoDS_Shape& shape = get(kernel, shapeId);
    BRepBuilderAPI_Copy copy(shape,true,false);
    BRepTools::Clean(copy.Shape());
    BRepBuilderAPI_GTransform maker(copy.Shape(), transform, true);
    if (!maker.IsDone() || maker.Shape().IsNull() || !BRepCheck_Analyzer(maker.Shape()).IsValid()) {
        fail(op, no_solution, "scaling produced invalid geometry");
    }
    history.collectComposed(copy, maker, shape);
    return history.finish(store(kernel, maker.Shape()));
}

} // namespace BRepExt
