// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepOffsetExt_Reconstruct.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>

namespace loft {
namespace {
gp_Pln outwardPlane(const TopoDS_Face &f) {
    BRepAdaptor_Surface surface(f);
    gp_Pnt at;
    gp_Vec du, dv;
    surface.D1((surface.FirstUParameter() + surface.LastUParameter()) / 2,
               (surface.FirstVParameter() + surface.LastVParameter()) / 2, at, du, dv);
    gp_Dir normal(du.Crossed(dv));
    if (f.Orientation() == TopAbs_REVERSED)
        normal.Reverse();
    return gp_Pln(at, normal);
}
TopoDS_Shape simplify(const TopoDS_Shape &shape) {
    ShapeUpgrade_UnifySameDomain unify(shape, true, true, false);
    unify.Build();
    return unify.Shape();
}
bool usable(const TopoDS_Shape &shape) {
    if (shape.IsNull() || !TopExp_Explorer(shape, TopAbs_SOLID).More() ||
        !BRepCheck_Analyzer(shape).IsValid())
        return false;
    return !intersectsItself(shape);
}
} // namespace
TopoDS_Shape reconstructOffset(const TopoDS_Shape &solid, const std::vector<TopoDS_Shape> &selected,
                               double distance, double tolerance, Lineage &faces) {
    if (selected.size() != 1)
        return {};
    const auto moved = TopoDS::Face(selected[0]);
    BRepAdaptor_Surface surface(moved);
    Adjacency adjacency;
    TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, adjacency);
    TopoDS_Shape result;
    // A full cylindrical cavity bounded only by perpendicular planes may break
    // through surrounding walls. Recut its exact finite analytic cylinder.
    if (surface.GetType() == GeomAbs_Cylinder && moved.Orientation() == TopAbs_REVERSED && distance < 0 &&
        std::abs(surface.LastUParameter() - surface.FirstUParameter() - 2 * M_PI) < loft::tolerance::angularResolution) {
        const auto c = surface.Cylinder();
        bool perpendicular = true;
        for (TopExp_Explorer e(moved, TopAbs_EDGE); e.More(); e.Next())
            for (const auto &f : adjacency.FindFromKey(e.Current()))
                if (!f.IsSame(moved)) {
                    BRepAdaptor_Surface s(TopoDS::Face(f));
                    if (s.GetType() != GeomAbs_Plane ||
                        std::abs(s.Plane().Axis().Direction().Dot(c.Axis().Direction())) < 1 - loft::tolerance::angle)
                        perpendicular = false;
                }
        if (perpendicular) {
            const auto start =
                c.Location().Translated(gp_Vec(c.Axis().Direction()) * surface.FirstVParameter());
            const auto tool =
                BRepPrimAPI_MakeCylinder(gp_Ax2(start, c.Axis().Direction()), c.Radius() - distance,
                                         surface.LastVParameter() - surface.FirstVParameter())
                    .Shape();
            Untouched<BRepAlgoAPI_Cut> cut(solid, tool);
            if (cut.IsDone())
                result = simplify(cut.Shape());
        }
    }
    if (surface.GetType() == GeomAbs_Plane) {
        const auto plane = outwardPlane(moved);
        bool convex = true;
        std::vector<gp_Pln> planes;
        for (TopExp_Explorer it(solid, TopAbs_FACE); it.More() && convex; it.Next()) {
            const auto f = TopoDS::Face(it.Current());
            if (BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane) {
                convex = false;
                break;
            }
            auto p = outwardPlane(f);
            for (TopExp_Explorer v(solid, TopAbs_VERTEX); v.More(); v.Next())
                if (gp_Vec(p.Location(), BRep_Tool::Pnt(TopoDS::Vertex(v.Current())))
                        .Dot(gp_Vec(p.Axis().Direction())) > tolerance * 10)
                    convex = false;
            if (f.IsSame(moved))
                p.Translate(gp_Vec(p.Axis().Direction()) * distance);
            planes.push_back(p);
        }
        if (convex) {
            Bnd_Box b;
            BRepBndLib::Add(solid, b);
            double x0, y0, z0, x1, y1, z1;
            b.Get(x0, y0, z0, x1, y1, z1);
            const double margin = std::abs(distance) + std::sqrt(b.SquareExtent()) + 1;
            result = BRepPrimAPI_MakeBox(gp_Pnt(x0 - margin, y0 - margin, z0 - margin),
                                         gp_Pnt(x1 + margin, y1 + margin, z1 + margin))
                         .Shape();
            for (const auto &p : planes) {
                auto boundary = BRepBuilderAPI_MakeFace(p).Face();
                auto half = BRepPrimAPI_MakeHalfSpace(boundary,
                                                      p.Location().Translated(-gp_Vec(p.Axis().Direction())))
                                .Solid();
                BRepAlgoAPI_Common common(result, half);
                if (!common.IsDone())
                    return {};
                result = common.Shape();
            }
            result = simplify(result);
        } else if (distance > 0) {
            // Extruding a planar face is exactly a face offset only when all
            // its adjacent supports are perpendicular planes. This allows a
            // concave step to disappear; it is not used across angled supports.
            bool perpendicular = true;
            for (TopExp_Explorer e(moved, TopAbs_EDGE); e.More(); e.Next())
                for (const auto &f : adjacency.FindFromKey(e.Current()))
                    if (!f.IsSame(moved)) {
                        BRepAdaptor_Surface s(TopoDS::Face(f));
                        if (s.GetType() != GeomAbs_Plane ||
                            std::abs(s.Plane().Axis().Direction().Dot(plane.Axis().Direction())) > loft::tolerance::angle)
                            perpendicular = false;
                    }
            if (perpendicular) {
                auto prism =
                    BRepPrimAPI_MakePrism(moved, gp_Vec(plane.Axis().Direction()) * distance).Shape();
                Untouched<BRepAlgoAPI_Fuse> fuse(solid, prism);
                if (fuse.IsDone())
                    result = simplify(fuse.Shape());
            }
        }
    }
    if (!usable(result))
        return {};
    // These reconstructions retain analytic plane/cylinder supports. Match
    // their loci, allowing one-to-many successors when a hole breaks out.
    for (TopExp_Explorer it(solid, TopAbs_FACE); it.More(); it.Next()) {
        const auto old = TopoDS::Face(it.Current());
        BRepAdaptor_Surface a(old);
        std::vector<TopoDS_Shape> successors;
        for (TopExp_Explorer jt(result, TopAbs_FACE); jt.More(); jt.Next()) {
            const auto next = TopoDS::Face(jt.Current());
            BRepAdaptor_Surface b(next);
            bool same = old.IsSame(next);
            if (a.GetType() == GeomAbs_Plane && b.GetType() == GeomAbs_Plane) {
                auto p = outwardPlane(old), q = outwardPlane(next);
                if (old.IsSame(moved))
                    p.Translate(gp_Vec(p.Axis().Direction()) * distance);
                same = p.Distance(q.Location()) <= tolerance &&
                       p.Axis().Direction().Dot(q.Axis().Direction()) > 1 - loft::tolerance::angle;
            } else if (a.GetType() == GeomAbs_Cylinder && b.GetType() == GeomAbs_Cylinder) {
                auto p = a.Cylinder(), q = b.Cylinder();
                same = std::abs(p.Radius() - (old.IsSame(moved) ? distance : 0) - q.Radius()) <= tolerance &&
                       std::abs(p.Axis().Direction().Dot(q.Axis().Direction())) > 1 - loft::tolerance::angle &&
                       gp_Lin(p.Axis()).Distance(q.Location()) <= tolerance;
            }
            if (same)
                successors.push_back(next);
        }
        faces.set(old, successors);
    }
    return solidIfSingle(result);
}
} // namespace loft
