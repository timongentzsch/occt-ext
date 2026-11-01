// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepOffsetExt_Parts.hxx>

#include <BRepSurfaceExt_FaceSeparation.hxx>

#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS_Compound.hxx>
#include <ChFi3d.hxx>
#include <IntTools_Context.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

#include <cmath>
#include <limits>

namespace loft {

namespace {

/// The (u, v) inside `face` on a grid of `perSide` squared, cell centres; along lines across it
/// where the grid misses it (a ring, a slot).
std::vector<gp_Pnt2d> parametersInside(const TopoDS_Face &face, int perSide) {
    const BRepAdaptor_Surface surface(face);
    BRepTopAdaptor_FClass2d inside(face, Precision::Confusion());
    const double u0 = surface.FirstUParameter(), u1 = surface.LastUParameter();
    const double v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
    std::vector<gp_Pnt2d> out;
    for (int i = 0; i < perSide; i++)
        for (int j = 0; j < perSide; j++) {
            const gp_Pnt2d uv(u0 + (u1 - u0) * (i + 0.5) / perSide, v0 + (v1 - v0) * (j + 0.5) / perSide);
            if (inside.Perform(uv) == TopAbs_IN)
                out.push_back(uv);
        }
    if (out.empty())
        out = parametersAcross(face, perSide, new IntTools_Context);
    return out;
}

NCollection_List<TopoDS_Shape> listOf(const std::vector<TopoDS_Shape> &faces) {
    NCollection_List<TopoDS_Shape> list;
    for (const TopoDS_Shape &face : faces)
        list.Append(face);
    return list;
}

} // namespace

OffsetFate classifyForOffset(const TopoDS_Shape &body, double distance) {
    OffsetFate fate;
    for (TopExp_Explorer faces(body, TopAbs_FACE); faces.More(); faces.Next()) {
        const TopoDS_Face &face = TopoDS::Face(faces.Current());
        const BRepAdaptor_Surface surface(face);
        if (surface.GetType() == GeomAbs_Plane)
            continue;
        // Along the surface's own normal the offset is this far.
        const double along = face.Orientation() == TopAbs_REVERSED ? -distance : distance;
        int collapsed = 0, measured = 0;
        // 9 x 9: a blend is the same all along, and a spline face that is tight in a part
        // smaller than a ninth of its width is left to the checks on the result.
        for (const gp_Pnt2d &uv : parametersInside(face, 9)) {
            BRepLProp_SLProps at(surface, uv.X(), uv.Y(), 2, Precision::Confusion());
            if (!at.IsCurvatureDefined())
                continue;
            measured++;
            // A principal curvature k becomes k / (1 - along * k): gone where along * k is 1.
            if (std::max(along * at.MaxCurvature(), along * at.MinCurvature()) >= 1.0 - collapseSlack)
                collapsed++;
        }
        if (collapsed == 0)
            continue;
        (collapsed == measured ? fate.vanish : fate.partial).push_back(face);
    }
    return fate;
}

bool intersectsItself(const TopoDS_Shape &shape) {
    // Only what the separation pre-pass cannot prove clear goes to the check (face_separation.h).
    for (const TopoDS_Shape &part : unprovenParts(shape)) {
        BOPAlgo_ArgumentAnalyzer analyzer;
        analyzer.SetShape1(part);
        analyzer.ArgumentTypeMode() = false;
        analyzer.SmallEdgeMode() = false;
        analyzer.SelfInterMode() = true;
        analyzer.Perform();
        if (analyzer.HasFaulty())
            return true;
    }
    return false;
}

std::vector<gp_Pnt> pointsInside(const TopoDS_Face &face, int perSide) {
    const BRepAdaptor_Surface surface(face);
    std::vector<gp_Pnt> out;
    for (const gp_Pnt2d &uv : parametersInside(face, perSide))
        out.push_back(surface.Value(uv.X(), uv.Y()));
    return out;
}

double clearance(const std::vector<TopoDS_Shape> &faces, const TopoDS_Shape &body,
                 const std::vector<TopoDS_Shape> &opened, double wall, bool wholeFaces) {
    // To its faces: a point inside a solid is at no distance from the solid itself.
    TopoDS_Compound skin, openings;
    BRep_Builder().MakeCompound(skin);
    BRep_Builder().MakeCompound(openings);
    for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next()) {
        bool open = false;
        for (const TopoDS_Shape &face : opened)
            open = open || face.IsSame(all.Current());
        BRep_Builder().Add(open ? openings : skin, all.Current());
    }
    double least = std::numeric_limits<double>::infinity();
    DistanceTo toSkin(skin), toOpening(openings);
    if (wholeFaces && opened.empty()) {
        // Measure whole trimmed faces, including their boundaries. Cell-centre samples miss
        // narrow clearance failures at a vanished concave blend or at an arc join.
        for (const TopoDS_Shape &face : faces) {
            try {
                const auto distance = toSkin.from(face);
                if (!distance)
                    return std::numeric_limits<double>::quiet_NaN();
                least = std::min(least, *distance);
            } catch (const Standard_Failure &) {
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
        return least;
    }
    for (const TopoDS_Shape &face : faces)
        for (const gp_Pnt &point : pointsInside(TopoDS::Face(face), 3)) {
            if (!opened.empty()) {
                // Beside an opening the wall ends, and a blend around the opening comes nearer.
                const auto beside = toOpening.from(point);
                if (beside && *beside < wall)
                    continue;
            }
            if (const auto distance = toSkin.from(point))
                least = std::min(least, *distance);
        }
    return least;
}

bool tangentAlong(const TopoDS_Edge &edge, const TopoDS_Face &a, const TopoDS_Face &b) {
    return !BRep_Tool::Degenerated(edge) && ChFi3d::IsTangentFaces(edge, a, b);
}

Lineage::Lineage(const TopoDS_Shape &input) {
    TopExp::MapShapes(input, TopAbs_FACE, input_);
    for (int i = 1; i <= input_.Extent(); i++)
        now_.push_back({input_.FindKey(i)});
}

void Lineage::throughHistory(const Handle(BRepTools_History) & history) {
    step([&](const TopoDS_Shape &face, std::vector<TopoDS_Shape> &out) {
        const NCollection_List<TopoDS_Shape> &next = history->Modified(face);
        if (next.IsEmpty() && !history->IsRemoved(face))
            out.push_back(face);
        for (const TopoDS_Shape &successor : next)
            if (successor.ShapeType() == TopAbs_FACE)
                out.push_back(successor);
    });
}

void Lineage::set(const TopoDS_Shape &face, std::vector<TopoDS_Shape> now) {
    const int index = input_.FindIndex(face);
    if (index != 0)
        now_[index - 1] = std::move(now);
}

void Lineage::then(const Lineage &next) {
    step([&](const TopoDS_Shape &face, std::vector<TopoDS_Shape> &out) {
        if (!next.input_.Contains(face)) {
            out.push_back(face);
            return;
        }
        for (const TopoDS_Shape &successor : next.of(face))
            out.push_back(successor);
    });
}

void Lineage::within(const TopoDS_Shape &result) {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> faces;
    TopExp::MapShapes(result, TopAbs_FACE, faces);
    step([&](const TopoDS_Shape &face, std::vector<TopoDS_Shape> &out) {
        const int index = faces.FindIndex(face);
        if (index == 0)
            return;
        for (const TopoDS_Shape &seen : out)
            if (seen.IsSame(face))
                return;
        out.push_back(faces.FindKey(index));
    });
}

const std::vector<TopoDS_Shape> &Lineage::of(const TopoDS_Shape &face) const {
    static const std::vector<TopoDS_Shape> none;
    const int index = input_.FindIndex(face);
    return index == 0 ? none : now_[index - 1];
}

std::vector<TopoDS_Shape> Lineage::all() const {
    std::vector<TopoDS_Shape> out;
    for (const auto &faces : now_)
        out.insert(out.end(), faces.begin(), faces.end());
    return out;
}

void Lineage::modified(FaceHistory &history) const {
    for (int i = 1; i <= input_.Extent(); i++) {
        const TopoDS_Shape &face = input_.FindKey(i);
        const auto &now = now_[i - 1];
        if (!history.tracks(face) || (now.size() == 1 && now[0].IsSame(face)))
            continue;
        if (now.empty())
            history.deleted(face);
        else
            history.modified(face, listOf(now));
    }
}

void Lineage::generated(FaceHistory &history) const {
    for (int i = 1; i <= input_.Extent(); i++)
        if (history.tracks(input_.FindKey(i)))
            history.generated(input_.FindKey(i), listOf(now_[i - 1]));
}

} // namespace loft
