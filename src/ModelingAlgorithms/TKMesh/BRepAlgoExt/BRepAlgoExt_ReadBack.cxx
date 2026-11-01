// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Geometry as a file gives it back (read_back.h).

#include <BRepAlgoExt_ReadBack.hxx>

#include <BRep_CurveRepresentation.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Conic.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_Conic.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Map.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Ax22d.hxx>
#include <gp_Lin.hxx>

namespace BRepExt {

Handle(Geom_Curve) basisOf(Handle(Geom_Curve) curve) {
    while (const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve)) curve = trimmed->BasisCurve();
    return curve;
}
Handle(Geom2d_Curve) basisOf(Handle(Geom2d_Curve) curve) {
    while (const auto trimmed = Handle(Geom2d_TrimmedCurve)::DownCast(curve)) curve = trimmed->BasisCurve();
    return curve;
}
Handle(Geom_Surface) basisOf(Handle(Geom_Surface) surface) {
    while (const auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface)) surface = trimmed->BasisSurface();
    return surface;
}

/// A frame as BRepTools and BinTools read it back: each direction's coordinates normalised once
/// (SetCoord), the frame rebuilt from the axis and x, and y reversed where the stored y says so.
static gp_Ax3 asRead(const gp_Ax3& frame) {
    const auto read = [](const gp_Dir& d) {
        gp_Dir out;
        out.SetCoord(d.X(), d.Y(), d.Z());
        return out;
    };
    const gp_Dir n = read(frame.Direction()), x = read(frame.XDirection()), y = read(frame.YDirection());
    gp_Ax3 again(frame.Location(), n, x);
    if (y.DotCross(n, x) < 0) again.YReverse();
    return again;
}

gp_Ax3 fixedFrame(gp_Ax3 frame) {
    // Reading back need not have a fixed point near a frame (it can alternate between two), so
    // the x direction, free about the axis, is turned by steps of 1e-13 until the frame reads
    // back as it is: about half of all frames are fixed at once, and over 100000 random frames
    // none needed more than a few dozen steps.
    const gp_Dir axis = readBack(frame.Direction()), x = frame.XDirection(), side = axis.Crossed(x);
    const bool direct = frame.Direct();
    for (int step = 0; step <= 256; ++step) {
        const double turn = (step % 2 ? 1 : -1) * ((step + 1) / 2) * 1e-13;
        gp_Ax3 again(frame.Location(), axis, gp_Dir(gp_Vec(x) * std::cos(turn) + gp_Vec(side) * std::sin(turn)));
        if (!direct) again.YReverse();
        const gp_Ax3 read = asRead(again);
        if (same(read.Direction(), again.Direction()) && same(read.XDirection(), again.XDirection()) &&
            same(read.YDirection(), again.YDirection()))
            return again;
    }
    return asRead(frame);
}

void settle(const Handle(Geom_Curve)& stored) {
    const auto curve = basisOf(stored);
    if (const auto line = Handle(Geom_Line)::DownCast(curve)) line->SetDirection(readBack(line->Lin().Direction()));
    else if (const auto conic = Handle(Geom_Conic)::DownCast(curve))
        for (int i = 0; i < 8; ++i) {
            const gp_Ax2 at = conic->Position(), again(at.Location(), readBack(at.Direction()), readBack(at.XDirection()));
            conic->SetPosition(again);
            if (same(again.Direction(), at.Direction()) && same(again.XDirection(), at.XDirection()) &&
                same(again.YDirection(), at.YDirection()))
                break;
        }
}

void settle(const Handle(Geom2d_Curve)& stored) {
    const auto curve = basisOf(stored);
    if (const auto line = Handle(Geom2d_Line)::DownCast(curve)) line->SetDirection(readBack(line->Direction()));
    else if (const auto conic = Handle(Geom2d_Conic)::DownCast(curve))
        for (int i = 0; i < 8; ++i) {
            const gp_Ax22d at = conic->Position(), again(at.Location(), readBack(at.XDirection()), readBack(at.YDirection()));
            conic->SetAxis(again);
            if (same(again.XDirection(), at.XDirection()) && same(again.YDirection(), at.YDirection())) break;
        }
}

void settleNew(const TopoDS_Shape& input, const TopoDS_Shape& result) {
    NCollection_Map<const Standard_Transient*> held;
    const auto curves = [](const TopoDS_Shape& edge, const auto& each) {
        const auto tedge = Handle(BRep_TEdge)::DownCast(edge.TShape());
        if (tedge.IsNull()) return;
        for (const auto& representation : tedge->Curves()) {
            if (representation->IsCurve3D() && !representation->Curve3D().IsNull()) each(representation->Curve3D());
            if (representation->IsCurveOnSurface()) each(representation->PCurve());
            if (representation->IsCurveOnClosedSurface()) each(representation->PCurve2());
        }
    };
    for (TopExp_Explorer faces(input, TopAbs_FACE); faces.More(); faces.Next()) {
        const auto surface = BRep_Tool::Surface(TopoDS::Face(faces.Current()));
        held.Add(surface.get());
        held.Add(basisOf(surface).get());
    }
    for (TopExp_Explorer edges(input, TopAbs_EDGE); edges.More(); edges.Next())
        curves(edges.Current(), [&](const auto& curve) {
            held.Add(curve.get());
            held.Add(basisOf(curve).get());
        });
    for (TopExp_Explorer faces(result, TopAbs_FACE); faces.More(); faces.Next()) {
        const auto surface = BRep_Tool::Surface(TopoDS::Face(faces.Current()));
        const auto elementary = Handle(Geom_ElementarySurface)::DownCast(basisOf(surface));
        if (elementary.IsNull() || held.Contains(surface.get()) || held.Contains(elementary.get())) continue;
        elementary->SetPosition(fixedFrame(elementary->Position()));
    }
    for (TopExp_Explorer edges(result, TopAbs_EDGE); edges.More(); edges.Next())
        curves(edges.Current(), [&](const auto& curve) {
            if (!held.Contains(curve.get()) && !held.Contains(basisOf(curve).get())) settle(curve);
        });
}

} // namespace BRepExt
