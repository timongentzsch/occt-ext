// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepAlgoExt_Geometry.hxx>

#include <BRepTools.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dHatch_Hatcher.hxx>
#include <Geom2d_Line.hxx>
#include <HatchGen_Domain.hxx>
#include <IntTools_Context.hxx>
#include <Standard_Failure.hxx>

namespace loft {

std::vector<gp_Pnt2d> parametersAcross(const TopoDS_Face& face, int lines, const Handle(IntTools_Context)& context) {
    std::vector<gp_Pnt2d> out;
    try {
        double u0, u1, v0, v1;
        BRepTools::UVBounds(face, u0, u1, v0, v1);
        if (Precision::IsInfinite(u0) || Precision::IsInfinite(u1) || Precision::IsInfinite(v0) ||
            Precision::IsInfinite(v1))
            return out;
        Geom2dHatch_Hatcher& hatcher = context->Hatcher(face);
        for (const bool alongV : {true, false})
            for (int i = 0; i < lines; ++i) {
                const double at = (i + 0.5) / lines;
                const Handle(Geom2d_Line) line =
                    alongV ? new Geom2d_Line(gp_Pnt2d(u0 + (u1 - u0) * at, v0), gp_Dir2d(0, 1))
                           : new Geom2d_Line(gp_Pnt2d(u0, v0 + (v1 - v0) * at), gp_Dir2d(1, 0));
                hatcher.ClrHatchings();
                const int hatching = hatcher.AddHatching(Geom2dAdaptor_Curve(line));
                hatcher.Trim();
                if (hatcher.TrimDone(hatching)) {
                    hatcher.ComputeDomains(hatching);
                    if (hatcher.IsDone(hatching))
                        for (int d = 1; d <= hatcher.NbDomains(hatching); ++d) {
                            const HatchGen_Domain& domain = hatcher.Domain(hatching, d);
                            if (!domain.HasFirstPoint() || !domain.HasSecondPoint()) continue;
                            const gp_Pnt2d uv =
                                line->Value((domain.FirstPoint().Parameter() + domain.SecondPoint().Parameter()) / 2);
                            if (context->IsPointInFace(face, uv)) out.push_back(uv);
                        }
                }
                hatcher.RemHatching(hatching);
            }
    } catch (const Standard_Failure&) {
    }
    return out;
}

} // namespace loft
