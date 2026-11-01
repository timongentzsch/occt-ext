// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomFill_Frenet.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>

namespace BRepExt {
std::vector<double> curveLocations(OcctKernel &kernel, uint32_t wireId, const std::vector<double> &fractions,
                                   bool transport) {
    const auto &shape = get(kernel, wireId);
    if (shape.ShapeType() != TopAbs_WIRE)
        fail("curveLocations", invalid_input, "path must be a wire");
    Handle(Adaptor3d_Curve) curve = new BRepAdaptor_CompCurve(TopoDS::Wire(shape), true);
    Handle(GeomFill_TrihedronLaw) frame;
    if (!transport) {
        frame = new GeomFill_Frenet();
        frame->SetCurve(curve);
    }
    const double length =
        GCPnts_AbscissaPoint::Length(*curve, curve->FirstParameter(), curve->LastParameter(), 1e-8);
    // Discrete parallel transport by shortest rotations. Unlike Frenet this has
    // no curvature-normal sign change at an inflection. Stations are independent
    // of requested array locations, so adding a copy cannot change earlier
    // frames.
    constexpr int steps = 1024;
    std::vector<gp_Vec> tangents, normals;
    auto tangentAt = [&](double u) {
        gp_Pnt p;
        gp_Vec t;
        curve->D1(u, p, t);
        if (t.Magnitude() < 1e-14)
            fail("curveLocations", invalid_input, "path has a singular tangent");
        return t.Normalized();
    };
    auto carry = [&](gp_Vec n, const gp_Vec &from, const gp_Vec &to) {
        const double dot = from.Dot(to);
        if (dot < -1 + 1e-10)
            fail("curveLocations", invalid_input, "rotation-minimising frame is undefined at a cusp");
        const auto axis = from.Crossed(to);
        n += axis.Crossed(n) + axis.Crossed(axis.Crossed(n)) / (1 + dot);
        return (n - to * n.Dot(to)).Normalized();
    };
    if (transport) {
        for (int i = 0; i <= steps; ++i) {
            const auto t = tangentAt(curve->FirstParameter() +
                                     (curve->LastParameter() - curve->FirstParameter()) * i / steps);
            gp_Vec n;
            if (i == 0) {
                gp_Vec reference = std::abs(t.Z()) < 0.9 ? gp_Vec(0, 0, 1) : gp_Vec(0, 1, 0);
                n = (reference - t * reference.Dot(t)).Normalized();
            } else
                n = carry(normals.back(), tangents.back(), t);
            tangents.push_back(t);
            normals.push_back(n);
        }
    }
    std::vector<double> values;
    for (double fraction : fractions) {
        if (!std::isfinite(fraction) || fraction < 0 || fraction > 1)
            fail("curveLocations", invalid_input, "fractions must be between zero and one");
        double u = curve->FirstParameter();
        if (fraction == 1)
            u = curve->LastParameter();
        else if (fraction > 0) {
            GCPnts_AbscissaPoint locate(1e-8, *curve, length * fraction, u);
            if (!locate.IsDone())
                fail("curveLocations", no_solution, "path length solve failed");
            u = locate.Parameter();
        }
        gp_Vec t, n, b;
        if (transport) {
            const int station =
                std::min(steps, int((u - curve->FirstParameter()) /
                                    (curve->LastParameter() - curve->FirstParameter()) * steps));
            t = tangentAt(u);
            n = carry(normals[station], tangents[station], t);
            b = t.Crossed(n);
        } else if (!frame->D0(u, t, n, b))
            fail("curveLocations", no_solution, "path frame is undefined");
        const auto p = curve->Value(u);
        values.insert(values.end(),
                      {p.X(), p.Y(), p.Z(), t.X(), t.Y(), t.Z(), n.X(), n.Y(), n.Z(), b.X(), b.Y(), b.Z()});
    }
    return values;
}
} // namespace BRepExt
