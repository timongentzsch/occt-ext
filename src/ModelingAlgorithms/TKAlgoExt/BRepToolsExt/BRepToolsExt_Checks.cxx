// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>

#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include <stdexcept>
#include <unordered_set>

namespace BRepExt {

} // namespace BRepExt
namespace BRepExt {

std::vector<TopoDS_Shape> subShapes(const OcctKernel& kernel, const char* operation,
                                    const TopoDS_Shape& owner, const std::vector<uint32_t>& ids,
                                    TopAbs_ShapeEnum type) {
    const char* noun = type == TopAbs_FACE ? "face" : type == TopAbs_VERTEX ? "vertex" : "edge";
    if (ids.empty()) fail(operation, invalid_input, std::string("no ") + noun + " given");

    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> owned;
    TopExp::MapShapes(owner, type, owned);

    std::vector<TopoDS_Shape> found;
    found.reserve(ids.size());
    std::unordered_set<int> seen;
    for (uint32_t id : ids) {
        const int index = owned.FindIndex(get(kernel, id));
        if (index == 0) {
            fail(operation, invalid_input,
                 std::string("selection is not ") + (type == TopAbs_FACE ? "a " : "an ") + noun +
                     " of the shape being modified");
        }
        if (!seen.insert(index).second) {
            fail(operation, invalid_input, std::string(noun) + " is given twice");
        }
        found.push_back(owned.FindKey(index));
    }
    return found;
}

std::string invalidity(const TopoDS_Shape& shape) {
    const BRepCheck_Analyzer check(shape);
    if (check.IsValid()) return {};
    std::string what;
    for (const auto kind : {TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SHELL, TopAbs_SOLID}) {
        NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> all;
        TopExp::MapShapes(shape, kind, all);
        for (int k = 1; k <= all.Extent(); ++k) {
            const auto& found = check.Result(all(k));
            if (found.IsNull()) continue;
            for (const auto status : found->Status())
                if (status != BRepCheck_NoError && what.size() < 300)
                    what += " " + std::to_string(int(kind)) + ":" + std::to_string(int(status));
            for (found->InitContextIterator(); found->MoreShapeInContext(); found->NextShapeInContext())
                for (const auto status : found->StatusOnShape())
                    if (status != BRepCheck_NoError && what.size() < 300)
                        what += " " + std::to_string(int(kind)) + "in" + std::to_string(int(found->ContextualShape().ShapeType())) + ":" +
                                std::to_string(int(status));
        }
    }
    return what.empty() ? " unknown" : what;
}

} // namespace BRepExt
