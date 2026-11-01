// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>

namespace loft {
// Keep every unselected dividing edge. Face dissolution permits only edges incident to
// a selected face; sharp edges still cannot be merged by same-domain unification.
EvolutionData cleanTopology(OcctKernel& k, uint32_t id, const std::vector<uint32_t>& selected,
                            const std::vector<int>& hashes, int bound) {
    const char* op="cleanTopology";
    const auto& shape=get(k,id);
    NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> edges,allowed,parts;
    TopExp::MapShapes(shape,TopAbs_EDGE,edges);TopExp::MapShapes(shape,parts);
    for(auto ref:selected) {
        const auto& part=get(k,ref);
        if(!parts.Contains(part) || (part.ShapeType()!=TopAbs_EDGE && part.ShapeType()!=TopAbs_FACE))
            fail(op,invalid_input,"selected edges or faces must belong to this body");
        if(part.ShapeType()==TopAbs_EDGE) allowed.Add(part);
        else TopExp::MapShapes(part,TopAbs_EDGE,allowed);
    }
    ShapeUpgrade_UnifySameDomain unify(shape,true,true,false);
    unify.SetSafeInputMode(true);
    if(!selected.empty()) for(const auto& edge:edges) if(!allowed.Contains(edge)) unify.KeepShape(edge);
    unify.Build();
    if(unify.Shape().IsNull() || !BRepCheck_Analyzer(unify.Shape()).IsValid())
        fail(op,no_solution,"redundant topology removal produced an invalid body");
    if(!selected.empty()) {
        NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> after;
        TopExp::MapShapes(unify.Shape(),TopAbs_EDGE,after);
        // A sharp/boundary edge cannot be dissolved without changing the geometry.
        for(auto ref:selected) {
            const auto& part=get(k,ref);
            if(part.ShapeType()==TopAbs_EDGE && after.Contains(part))
                fail(op,unsupported_geometry,"selected edge is a boundary or separates different supports");
        }
        if(after.Extent()==edges.Extent()) fail(op,unsupported_geometry,"selected faces have no redundant dividing edges");
    }
    FaceHistory history(op,hashes,bound);auto changes=unify.History();
    if(!changes.IsNull()) for(TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next()) {
        const auto& face=it.Current();if(!history.tracks(face)) continue;
        history.modified(face,changes->Modified(face));
        history.generated(face,changes->Generated(face));
        if(changes->IsRemoved(face)) history.deleted(face);
    }
    return history.finish(store(k,unify.Shape()));
}
}
