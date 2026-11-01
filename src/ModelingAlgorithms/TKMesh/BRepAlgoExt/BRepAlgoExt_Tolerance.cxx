// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepAlgoExt_Tolerance.hxx>

#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <cmath>
#include <limits>

namespace loft::tolerance {

double own(const TopoDS_Shape& shape) {
    switch (shape.ShapeType()) {
    case TopAbs_VERTEX: return BRep_Tool::Tolerance(TopoDS::Vertex(shape));
    case TopAbs_EDGE: return BRep_Tool::Tolerance(TopoDS::Edge(shape));
    case TopAbs_FACE: return BRep_Tool::Tolerance(TopoDS::Face(shape));
    default: return 0;
    }
}

double largest(const TopoDS_Shape& shape) {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> seen;
    TopExp::MapShapes(shape, seen);
    double value = 0;
    for (int i = 1; i <= seen.Extent(); ++i) {
        const double local = own(seen(i));
        if (!std::isfinite(local)) return std::numeric_limits<double>::infinity();
        value = std::max(value, local);
    }
    return value;
}

} // namespace loft::tolerance
