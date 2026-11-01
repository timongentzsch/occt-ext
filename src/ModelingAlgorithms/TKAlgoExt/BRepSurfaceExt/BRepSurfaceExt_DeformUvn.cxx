// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepSurfaceExt_DeformUvn.hxx>
#include <GeomAlgoExt_CubicThrough.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepLib_FindSurface.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <NCollection_Array2.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax3.hxx>
#include <gp_Lin.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <algorithm>
#include <cmath>

namespace loft {
DeformUVN::DeformUVN(const char *operation, const TopoDS_Face &reference, const TopoDS_Face &target,
                     Vec3 s, Vec3 o, bool sw, bool mi, bool normal)
    : op(operation), surface(BRep_Tool::Surface(target)), scale(s), offset(o), swap(sw), mirror(mi) {
    BRepAdaptor_Surface flat(reference);
    if (flat.GetType() != GeomAbs_Plane)
        fail(op, unsupported_geometry, "reference face must be planar");
    plane = flat.Plane();
    for (double value : {s.x, s.y, s.z, o.x, o.y, o.z})
        if (!std::isfinite(value)) fail(op, invalid_input, "scale and offset must be finite");
    if (s.x == 0 || s.y == 0) fail(op, invalid_input, "UV scale must be nonzero");
    BRepTools::UVBounds(reference, a, b, c, d);
    BRepTools::UVBounds(target, u0, u1, v0, v1);
    for (double value : {a, b, c, d, u0, u1, v0, v1})
        if (!std::isfinite(value)) fail(op, invalid_input, "finite face parameter domains required");
    if (b <= a || d <= c || u1 <= u0 || v1 <= v0)
        fail(op, invalid_input, "nondegenerate face parameter domains required");
    up = (plane.Direct() ? 1.0 : -1.0) * (reference.Orientation() == TopAbs_REVERSED ? -1.0 : 1.0);
    out = (target.Orientation() == TopAbs_REVERSED ? -1.0 : 1.0) * (normal ? -1.0 : 1.0);
}
gp_Pnt DeformUVN::local(const gp_Pnt &p) const {
    const gp_Vec r(plane.Location(), p);
    return {r.Dot(gp_Vec(plane.XAxis().Direction())), r.Dot(gp_Vec(plane.YAxis().Direction())),
            r.Dot(gp_Vec(plane.Axis().Direction())) * up};
}
gp_Pnt DeformUVN::world(const gp_Pnt &p) const {
    return plane.Location().Translated(gp_Vec(plane.XAxis().Direction()) * p.X() +
                                      gp_Vec(plane.YAxis().Direction()) * p.Y() +
                                      gp_Vec(plane.Axis().Direction()) * (p.Z() * up));
}
gp_XY DeformUVN::uv(const gp_Pnt &p) const {
    auto r = local(p);
    double u = (r.X() - a) / (b - a), v = (r.Y() - c) / (d - c);
    if (swap) std::swap(u, v);
    if (mirror) u = 1 - u;
    return {u0 + (u * scale.x + offset.x) * (u1 - u0),
            v0 + (v * scale.y + offset.y) * (v1 - v0)};
}
gp_Pnt DeformUVN::point(const gp_Pnt &p) const {
    auto at = uv(p);
    gp_Pnt on;
    gp_Vec du, dv;
    surface->D1(at.X(), at.Y(), on, du, dv);
    auto n = du.Crossed(dv);
    const double height = local(p).Z() * scale.z + offset.z;
    if (n.Magnitude() < 1e-14) {
        if (height != 0) fail(op, no_solution, "target normal is singular");
        return on;
    }
    auto result = on.Translated(n.Normalized() * (height * out));
    for (double v : {result.X(), result.Y(), result.Z()})
        if (!std::isfinite(v)) fail(op, no_solution, "UVN map is not finite");
    return result;
}
void DeformUVN::checkChart(double umin, double umax, double vmin, double vmax) const {
    // Periodic support surfaces are allowed only inside one explicit, open chart. Whole
    // cylinders and any crossing of the selected chart's cut require splitting/repicking.
    if (surface->IsUPeriodic() &&
        (u1 - u0 >= surface->UPeriod() - 1e-9 || umin < u0 - 1e-9 || umax > u1 + 1e-9))
        fail(op, unsupported_geometry, "periodic U seam crossing: split the target and moving body into open charts");
    if (surface->IsVPeriodic() &&
        (v1 - v0 >= surface->VPeriod() - 1e-9 || vmin < v0 - 1e-9 || vmax > v1 + 1e-9))
        fail(op, unsupported_geometry, "periodic V seam crossing: split the target and moving body into open charts");
    if (surface->Continuity() < GeomAbs_C2)
        fail(op, unsupported_geometry, "target continuity breaks require separate charts");
}
} // namespace loft

namespace loft {
// The faces that carry a deformation from one curve to another, for curves taken as the
// reference and the target: a flat strip along a
// straight reference, wide enough for what moves, and the same strip laid along the target by
// arc length, straight across from it in the target's plane. Deform then runs as between faces, so
// shares of the reference go to shares of the target and distances off it stay.
EvolutionData deformSupports(OcctKernel &kernel, uint32_t referenceId, uint32_t targetId, uint32_t movingId,
                             double tolerance) {
    const char *op = "deformSupports";
    if (!std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "tolerance must be positive");
    const auto wireOf = [&](uint32_t id) {
        const auto &shape = get(kernel, id);
        if (shape.ShapeType() == TopAbs_WIRE)
            return TopoDS::Wire(shape);
        if (shape.ShapeType() != TopAbs_EDGE)
            fail(op, invalid_input, "reference and target have to be curves");
        return BRepBuilderAPI_MakeWire(TopoDS::Edge(shape)).Wire();
    };
    BRepAdaptor_CompCurve reference(wireOf(referenceId), true), target(wireOf(targetId), true);
    const double referenceLength = GCPnts_AbscissaPoint::Length(reference);
    const double targetLength = GCPnts_AbscissaPoint::Length(target, tolerance / 1000);
    const gp_Pnt start = reference.Value(reference.FirstParameter()), end = reference.Value(reference.LastParameter());
    if (referenceLength <= tolerance || targetLength <= tolerance || start.Distance(end) <= tolerance)
        fail(op, invalid_input, "the reference and target have to be open curves with length");
    for (int i = 1; i < 16; ++i) {
        const double u = reference.FirstParameter() + (reference.LastParameter() - reference.FirstParameter()) * i / 16;
        if (gp_Lin(start, gp_Vec(start, end)).Distance(reference.Value(u)) > tolerance)
            fail(op, unsupported_geometry, "the reference curve has to be straight");
    }
    const gp_Dir running(gp_Vec(start, end));
    // The plane: what moves lies in it with the reference. Its side of the reference is the positive one.
    Bnd_Box box;
    BRepBndLib::Add(get(kernel, movingId), box);
    if (box.IsVoid())
        fail(op, invalid_input, "there is nothing to move");
    double x0, y0, z0, x1, y1, z1;
    box.Get(x0, y0, z0, x1, y1, z1);
    const gp_Pnt middle((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
    gp_Vec beside(start, middle);
    beside -= gp_Vec(running) * beside.Dot(gp_Vec(running));
    // The reference's plane: the plane of what moves, where it is flat and holds the reference's
    // direction; otherwise the target's, or else the one through the middle of what moves.
    BRepLib_FindSurface targetPlane(wireOf(targetId), tolerance, true);
    const auto square = [&](const gp_Vec &n) { return n - gp_Vec(running) * n.Dot(gp_Vec(running)); };
    gp_Vec normal;
    BRepLib_FindSurface movingPlane(get(kernel, movingId), tolerance, true);
    if (movingPlane.Found()) {
        const gp_Vec n(occ::down_cast<Geom_Plane>(movingPlane.Surface())->Pln().Axis().Direction());
        if (std::abs(n.Dot(gp_Vec(running))) <= 1e-9)
            normal = n;
    }
    if (normal.Magnitude() <= gp::Resolution() && targetPlane.Found())
        normal = square(gp_Vec(occ::down_cast<Geom_Plane>(targetPlane.Surface())->Pln().Axis().Direction()));
    if (normal.Magnitude() <= gp::Resolution())
        normal = gp_Vec(running).Crossed(beside);
    if (normal.Magnitude() <= tolerance)
        fail(op, invalid_input, "the reference and what moves do not fix a plane");
    normal.Normalize();
    // What moves lies on the positive side of the reference.
    if (gp_Vec(running).Crossed(beside).Dot(normal) < 0)
        normal.Reverse();
    const gp_Vec across = normal.Crossed(gp_Vec(running));
    // How far from the reference it reaches, across, and past which ends it may not go.
    double reach = 0;
    for (const double x : {x0, x1})
        for (const double y : {y0, y1})
            for (const double z : {z0, z1})
                reach = std::max(reach, std::abs(gp_Vec(start, gp_Pnt(x, y, z)).Dot(across)));
    reach = std::max(reach * 1.01, referenceLength * 1e-3);
    gp_Ax3 frame(start, gp_Dir(normal), running);
    const TopoDS_Face flat = BRepBuilderAPI_MakeFace(gp_Pln(frame), 0, referenceLength, -reach, reach).Face();
    // The target's plane, and from it the way across at each place. A straight target lies in the
    // reference's plane turned about it, or takes the normal as it stands.
    gp_Vec targetNormal = normal;
    if (targetPlane.Found()) {
        targetNormal = gp_Vec(occ::down_cast<Geom_Plane>(targetPlane.Surface())->Pln().Axis().Direction());
        if (targetNormal.Dot(normal) < 0)
            targetNormal.Reverse();
    } else {
        for (int i = 1; i < 16; ++i) {
            const gp_Pnt p0 = target.Value(target.FirstParameter()), p1 = target.Value(target.LastParameter());
            const double u = target.FirstParameter() + (target.LastParameter() - target.FirstParameter()) * i / 16;
            if (gp_Lin(p0, gp_Vec(p0, p1)).Distance(target.Value(u)) > tolerance)
                fail(op, unsupported_geometry, "the target curve has to lie in a plane");
        }
    }
    // The target is taken the way the reference runs, whichever way it was drawn.
    const bool against = gp_Vec(target.Value(target.FirstParameter()), target.Value(target.LastParameter())).Dot(gp_Vec(running)) < 0;
    const auto placed = [&](double share) {
        GCPnts_AbscissaPoint place(target, (against ? 1 - share : share) * targetLength, target.FirstParameter(), tolerance / 1000);
        if (!place.IsDone())
            fail(op, no_solution, "the target curve could not be walked along");
        gp_Pnt p;
        gp_Vec d;
        target.D1(place.Parameter(), p, d);
        gp_Vec side = targetNormal.Crossed(against ? d.Reversed() : d);
        if (side.Magnitude() <= gp::Resolution())
            fail(op, no_solution, "the target curve stands square to the plane");
        side.Normalize();
        return std::pair{p.Translated(side * -reach), p.Translated(side * reach)};
    };
    for (int n = 16;; n *= 2) {
        if (n > 4096)
            fail(op, no_solution, "the target curve bends too sharply to lay a strip along");
        std::vector<double> places;
        std::vector<gp_Pnt> low, high;
        for (int k = 0; k <= n; ++k) {
            places.push_back(double(k) / n);
            const auto [a, b] = placed(places.back());
            low.push_back(a);
            high.push_back(b);
        }
        const GeomAlgoExt::CubicThrough law(places, GeomAlgoExt::CubicThrough::Ends::Follow);
        const auto lowPoles = law.through(low), highPoles = law.through(high);
        NCollection_Array2<gp_Pnt> poles(1, law.poles, 1, 2);
        for (int k = 1; k <= law.poles; ++k) {
            poles(k, 1) = lowPoles[size_t(k - 1)];
            poles(k, 2) = highPoles[size_t(k - 1)];
        }
        NCollection_Array1<double> vKnots(1, 2);
        NCollection_Array1<int> vMults(1, 2);
        vKnots(1) = 0;
        vKnots(2) = 1;
        vMults(1) = vMults(2) = 2;
        Handle(Geom_BSplineSurface) surface = new Geom_BSplineSurface(poles, law.knots, vKnots, law.mults, vMults, 3, 1);
        double worst = 0;
        for (int k = 0; k < n; ++k) {
            const double u = (k + 0.5) / n;
            const auto [a, b] = placed(u);
            worst = std::max({worst, surface->Value(u, 0).Distance(a), surface->Value(u, 1).Distance(b)});
        }
        if (worst > tolerance)
            continue;
        TopoDS_Compound both;
        BRep_Builder builder;
        builder.MakeCompound(both);
        builder.Add(both, flat);
        builder.Add(both, BRepBuilderAPI_MakeFace(surface, Precision::Confusion()).Face());
        return stored(kernel, both);
    }
}
} // namespace loft
