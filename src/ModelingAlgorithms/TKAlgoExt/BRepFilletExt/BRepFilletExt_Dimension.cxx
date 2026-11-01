// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAlgoExt_RoundFace.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <cmath>
#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <gp_Pln.hxx>

namespace BRepExt {
namespace {
bool sameSupport(const TopoDS_Face& a, const TopoDS_Face& b, double tolerance) {
    BRepAdaptor_Surface surface(a);
    if (surface.GetType() != BRepAdaptor_Surface(b).GetType()) return false;
    const auto other = BRep_Tool::Surface(b);
    if (BRep_Tool::Surface(a) == other) return true;
    for (double f : {.2, .5, .8}) for (double g : {.2, .5, .8}) {
        const auto p = surface.Value(surface.FirstUParameter() + f * (surface.LastUParameter() - surface.FirstUParameter()),
                                     surface.FirstVParameter() + g * (surface.LastVParameter() - surface.FirstVParameter()));
        GeomAPI_ProjectPointOnSurf projection(p, other);
        if (!projection.NbPoints() || projection.LowerDistance() > tolerance) return false;
        double u, v;projection.LowerDistanceParameters(u,v);
        gp_Pnt at;gp_Vec a0,a1,b0,b1;
        surface.D1(surface.FirstUParameter() + f * (surface.LastUParameter() - surface.FirstUParameter()),
                   surface.FirstVParameter() + g * (surface.LastVParameter() - surface.FirstVParameter()),at,a0,a1);
        other->D1(u,v,at,b0,b1);
        auto na=a0.Crossed(a1),nb=b0.Crossed(b1);
        if (na.Magnitude()<BRepExt::tolerance::length || nb.Magnitude()<BRepExt::tolerance::length ||
            std::abs(na.Normalized().Dot(nb.Normalized()))<1-BRepExt::tolerance::angle) return false;
    }
    return true;
}
} // namespace
RoundFace recognizeRound(const TopoDS_Shape &shape, const TopoDS_Face &face, double tolerance) {
    RoundFace b;
    BRepAdaptor_Surface s(face);
    if (s.GetType() != GeomAbs_Cylinder && s.GetType() != GeomAbs_Torus)
        return b;
    if (s.GetType() == GeomAbs_Cylinder && s.LastUParameter() - s.FirstUParameter() >= 2 * M_PI - BRepExt::tolerance::angle) return b;
    b.toroidal = s.GetType() == GeomAbs_Torus;
    if (b.toroidal) b.torus = s.Torus();
    else b.cylinder = s.Cylinder();
    Adjacency adjacency;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacency);
    ShapeSet neighbours;
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
        for (const auto &other : adjacency.FindFromKey(e.Current()))
            if (!other.IsSame(face))
                neighbours.Add(other);
    for (const auto &neighbour : neighbours) {
        bool tangent = false, bad = false;
        for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next()) {
            bool shared = false;
            for (TopExp_Explorer other(neighbour, TopAbs_EDGE); other.More(); other.Next())
                shared |= e.Current().IsSame(other.Current());
            if (!shared || BRep_Tool::Degenerated(TopoDS::Edge(e.Current()))) continue;
            const auto edge = TopoDS::Edge(e.Current());
            BRepAdaptor_Curve curve(edge);
            bool smooth = true;
            for (int i = 0; i < 17; ++i) {
                const double t = curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * (i + .5) / 17;
                const auto a = surfaceJet("filletDimension", edge, face, t);
                const auto c = surfaceJet("filletDimension", edge, TopoDS::Face(neighbour), t);
                smooth &= a.p.Distance(c.p) <= tolerance && a.n.Dot(c.n) > 1 - BRepExt::tolerance::angle;
            }
            tangent |= smooth;
            bad |= !smooth;
        }
        if (tangent && !bad) b.supports.push_back(TopoDS::Face(neighbour));
    }
    b.valid = b.supports.size() == 2;
    return b;
}
std::vector<double> filletDimension(OcctKernel &kernel, uint32_t shapeId, uint32_t faceId, double tolerance) {
    const auto &shape = get(kernel, shapeId);
    const auto owned = subShapes(kernel, "filletDimension", shape, {faceId}, TopAbs_FACE);
    const auto face = TopoDS::Face(owned[0]);
    const auto b = recognizeRound(shape, face, tolerance);
    if (!b.valid)
        return {};
    BRepAdaptor_Surface s(face);
    const auto at = s.Value((s.FirstUParameter() + s.LastUParameter()) / 2,
                            (s.FirstVParameter() + s.LastVParameter()) / 2);
    gp_Pnt center;
    double radius;
    if (b.toroidal) {
        const auto axis = b.torus.Axis();
        const gp_Vec radial = gp_Vec(b.torus.Location(), at) - gp_Vec(axis.Direction()) * gp_Vec(b.torus.Location(), at).Dot(gp_Vec(axis.Direction()));
        if (radial.Magnitude() <= tolerance) return {};
        center = b.torus.Location().Translated(radial.Normalized() * b.torus.MajorRadius());
        radius = b.torus.MinorRadius();
    } else {
        const auto axis = b.cylinder.Axis();
        center = axis.Location().Translated(gp_Vec(axis.Direction()) * gp_Vec(axis.Location(), at).Dot(gp_Vec(axis.Direction())));
        radius = b.cylinder.Radius();
    }
    return {radius, center.X(), center.Y(), center.Z(), at.X(), at.Y(), at.Z()};
}
EvolutionData refillet(OcctKernel &kernel, uint32_t shapeId, uint32_t faceId, double radius, double tolerance,
                       const std::vector<int> &hashes, int bound) {
    const char *op = "refillet";
    const auto &shape = get(kernel, shapeId);
    const auto owned = subShapes(kernel, op, shape, {faceId}, TopAbs_FACE);
    const auto face = TopoDS::Face(owned[0]);
    const auto b = recognizeRound(shape, face, tolerance);
    if (!b.valid || !std::isfinite(radius) || radius <= tolerance)
        fail(op, invalid_input,
             "face must be a cylindrical or toroidal fillet tangent to two supports; "
             "radius must be positive");
    BRepAlgoAPI_Defeaturing remove;
    remove.SetShape(shape);
    remove.AddFaceToRemove(face);
    remove.Build(progress(true)->Start());
    checkCancelled(op);
    if (!remove.IsDone())
        fail(op, no_solution, "could not remove the existing blend");
    const auto sharp = remove.Shape();
    Adjacency adjacency;
    TopExp::MapShapesAndAncestors(sharp, TopAbs_EDGE, TopAbs_FACE, adjacency);
    BRepFilletAPI_MakeFillet maker(sharp);
    int count = 0;
    for (int i = 1; i <= adjacency.Extent(); ++i) {
        bool found[2] = {false, false};
        for (const auto &adjacent : adjacency.FindFromIndex(i)) {
            for (int j = 0; j < 2; ++j)
                if (sameSupport(b.supports[j], TopoDS::Face(adjacent), tolerance)) found[j] = true;
        }
        if (found[0] && found[1]) {
            maker.Add(radius, TopoDS::Edge(adjacency.FindKey(i)));
            ++count;
        }
    }
    if (count != 1)
        fail(op, no_solution, "removing this blend did not expose exactly one support edge");
    maker.Build();
    if (!maker.IsDone() || !BRepCheck_Analyzer(maker.Shape()).IsValid())
        fail(op, no_solution, "replacement fillet failed");
    FaceHistory history(op, hashes, bound);
    history.collectComposed(remove, maker, shape, [&](const TopoDS_Shape &original, NCollection_List<TopoDS_Shape> &final) {
        if (!original.IsSame(face))
            return;
        // The consumed fillet's successor is the blend generated by its recovered edge.
        for (int i = 1; i <= adjacency.Extent(); ++i)
            for (const auto &generated : maker.Generated(adjacency.FindKey(i)))
                if (generated.ShapeType() == TopAbs_FACE)
                    final.Append(generated);
    });
    return history.finish(store(kernel, maker.Shape()));
}
} // namespace BRepExt
