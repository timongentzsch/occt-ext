// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// A planar face beside a face of a cylinder, a cone, a sphere or a torus whose axis is normal to
// its plane: the wall of a hole or a boss on the face it passes through, a countersink, a rim
// blend. Such a plane cuts the surface in circles about the axis at its height h, in closed
// form, and the surface's points within a gap of the plane are those whose height is within that
// gap of h: a band about each circle no wider than gap / sin t, where t is the angle the two
// surfaces make along it, the same all round. So the two faces meet only along their shared
// edges, and come within their tolerances of each other only beside them, where every circle
//   - is a shared edge, whole, with the surfaces at least kCrease apart along it (the crease that
//     face_separation's NEAR_EDGE admits), or
//   - has its band clear of the box of either face.
// Anything else (an axis not normal to the plane, the apex of a cone, a plane through a sphere's
// pole or touching a torus, a circle that is part of an edge only) is left to the check.

#include <BRepSurfaceExt_FaceSeparationAnalytic.hxx>

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <BndLib.hxx>
#include <Bnd_Box.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Circ.hxx>

#include <cmath>
#include <optional>
#include <vector>

namespace loft {
namespace {

/// sin 15 degrees: two surfaces crossing at less than this come within a gap of each other
/// more than four gaps away from their line, more than face_separation's NEAR_EDGE allows.
constexpr double kCrease = 0.2588;

struct Ring {
    gp_Circ circle;
    double crease = 0; // the sine of the angle between the surfaces along it
};

bool revolution(const BRepAdaptor_Surface& surface) {
    const auto type = surface.GetType();
    return type == GeomAbs_Cylinder || type == GeomAbs_Cone || type == GeomAbs_Sphere || type == GeomAbs_Torus;
}

/// The axis of a surface of revolution, along `normal` for a sphere.
gp_Ax1 axisOf(const BRepAdaptor_Surface& surface, const gp_Dir& normal) {
    switch (surface.GetType()) {
    case GeomAbs_Cylinder: return surface.Cylinder().Axis();
    case GeomAbs_Cone: return surface.Cone().Axis();
    case GeomAbs_Torus: return surface.Torus().Axis();
    default: return gp_Ax1(surface.Sphere().Location(), normal);
    }
}

/// Where `plane` cuts the surface of `face`, about its axis normal to the plane; none where the
/// pair is not of that kind or the cut is singular.
std::optional<std::vector<Ring>> ringsOf(const gp_Pln& plane, const BRepAdaptor_Surface& surface, double gap) {
    const gp_Dir normal = plane.Axis().Direction();
    if (!revolution(surface)) return std::nullopt;
    const gp_Ax1 axis = axisOf(surface, normal);
    if (!axis.Direction().IsParallel(normal, Precision::Angular())) return std::nullopt;
    const gp_Dir up = axis.Direction();
    const double h = gp_Vec(axis.Location(), plane.Location()).Dot(gp_Vec(up));
    const gp_Ax2 at(axis.Location().Translated(gp_Vec(up) * h), up);
    std::vector<Ring> rings;
    const auto ring = [&](double radius, double crease) {
        if (!(radius > gap)) return false;
        rings.push_back({gp_Circ(at, radius), crease});
        return true;
    };
    switch (surface.GetType()) {
    case GeomAbs_Cylinder:
        if (!ring(surface.Cylinder().Radius(), 1)) return std::nullopt;
        break;
    case GeomAbs_Cone: {
        // At height h along the axis from its reference circle, R + h tan(a).
        const auto cone = surface.Cone();
        if (!ring(cone.RefRadius() + h * std::tan(cone.SemiAngle()), std::abs(std::cos(cone.SemiAngle()))))
            return std::nullopt;
        break;
    }
    case GeomAbs_Sphere: {
        const double r = surface.Sphere().Radius(), across = r * r - h * h;
        if (across <= 0 || !ring(std::sqrt(across), std::sqrt(across) / r)) return std::nullopt;
        break;
    }
    case GeomAbs_Torus: {
        const auto torus = surface.Torus();
        const double r = torus.MinorRadius(), across = r * r - h * h;
        if (across <= 0) return std::nullopt;
        const double s = std::sqrt(across);
        if (!ring(torus.MajorRadius() + s, s / r) || !ring(torus.MajorRadius() - s, s / r)) return std::nullopt;
        break;
    }
    default: return std::nullopt;
    }
    return rings;
}

/// Whether `edge` is the whole of `circle`, within its own tolerance: the two surfaces cross on
/// the circle, so there they cross on the edge.
bool wholeOf(const TopoDS_Edge& edge, const gp_Circ& circle) {
    const double gap = std::max(Precision::Confusion(), BRep_Tool::Tolerance(edge));
    TopoDS_Vertex first, last;
    TopExp::Vertices(edge, first, last);
    if (first.IsNull() || !first.IsSame(last)) return false;
    const BRepAdaptor_Curve curve(edge);
    if (curve.GetType() != GeomAbs_Circle) return false;
    const gp_Circ own = curve.Circle();
    return std::abs(curve.LastParameter() - curve.FirstParameter() - 2 * M_PI) <= Precision::PConfusion() &&
           own.Location().Distance(circle.Location()) <= gap && std::abs(own.Radius() - circle.Radius()) <= gap &&
           own.Axis().Direction().IsParallel(circle.Axis().Direction(), Precision::Angular());
}

bool planeBesideRevolution(const TopoDS_Face& flat, const TopoDS_Face& round,
                           const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& shared) {
    const double gap = tolerance::largest(flat) + tolerance::largest(round);
    const BRepAdaptor_Surface plane(flat), surface(round);
    const auto rings = ringsOf(plane.Plane(), surface, gap);
    if (!rings) return false;
    // Every shared edge a whole circle of the cut, every shared vertex one of such an edge.
    std::vector<char> matched(rings->size(), 0);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> ends;
    for (int k = 1; k <= shared.Extent(); ++k) {
        if (shared(k).ShapeType() != TopAbs_EDGE) continue;
        const auto edge = TopoDS::Edge(shared(k));
        bool found = false;
        for (size_t i = 0; i < rings->size() && !found; ++i)
            if (wholeOf(edge, (*rings)[i].circle)) found = matched[i] = 1;
        if (!found) return false;
        TopExp::MapShapes(edge, TopAbs_VERTEX, ends);
    }
    for (int k = 1; k <= shared.Extent(); ++k)
        if (shared(k).ShapeType() == TopAbs_VERTEX && !ends.Contains(shared(k))) return false;
    Bnd_Box flatBox, roundBox;
    BRepBndLib::Add(flat, flatBox, false);
    BRepBndLib::Add(round, roundBox, false);
    flatBox.Enlarge(gap);
    roundBox.Enlarge(gap);
    for (size_t i = 0; i < rings->size(); ++i) {
        const Ring& ring = (*rings)[i];
        if (ring.crease < kCrease) return false;
        if (matched[i]) continue;
        // A circle of the cut that is no edge of theirs: its band clear of one of the faces.
        Bnd_Box band;
        BndLib::Add(ring.circle, 0, band);
        band.Enlarge(gap / ring.crease + gap);
        if (!band.IsOut(flatBox) && !band.IsOut(roundBox)) return false;
    }
    return true;
}

} // namespace

bool analyticBeside(const TopoDS_Face& a, const TopoDS_Face& b,
                    const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& shared) {
    const BRepAdaptor_Surface sa(a, false), sb(b, false);
    if (sa.GetType() == GeomAbs_Plane && revolution(sb)) return planeBesideRevolution(a, b, shared);
    if (sb.GetType() == GeomAbs_Plane && revolution(sa)) return planeBesideRevolution(b, a, shared);
    return false;
}

bool analyticWorth(const TopoDS_Shape& shape) {
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> owners;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, owners);
    // An edge between a plane and a surface of revolution about an axis normal to it.
    for (int e = 1; e <= owners.Extent(); ++e) {
        std::optional<gp_Dir> normal;
        std::vector<BRepAdaptor_Surface> round;
        for (const auto& face : owners(e)) {
            const BRepAdaptor_Surface surface(TopoDS::Face(face), false);
            if (surface.GetType() == GeomAbs_Plane) normal = surface.Plane().Axis().Direction();
            else if (revolution(surface)) round.push_back(surface);
        }
        for (const auto& surface : round)
            if (normal && axisOf(surface, *normal).Direction().IsParallel(*normal, Precision::Angular())) return true;
    }
    return false;
}

} // namespace loft
