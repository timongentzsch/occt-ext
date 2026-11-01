// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

namespace BRepExt {
namespace surfacing {
struct Boundary {
    TopoDS_Edge edge;
    TopoDS_Face face;
    occ::handle<Geom_Curve> curve;
    double first, last;
    gp_Pnt point(double f) const;
    SurfaceJet jet(const char *op, double f) const;
};
struct Constraint {
    Boundary reference;
    gp_XY start, end;
    int order;
    // How far the surface's way off the boundary leans along it at its two ends, as the tangent of
    // the angle from straight across: towards the other end when positive. Between them it runs
    // evenly. Only where `skewed`. Unsourced convention: the origin and sense of the angle are chosen,
    // not derived.
    bool skewed = false;
    double skewStart = 0, skewEnd = 0;
    double skewAt(double f) const { return skewStart * (1 - f) - skewEnd * f; }
};
TopoDS_Face faceOf(OcctKernel &, const char *, uint32_t);
Boundary boundaryOf(OcctKernel &, const char *, uint32_t edge, uint32_t face);
occ::handle<Geom_BSplineSurface> editable(const char *, const TopoDS_Face &, bool fitting = true);
SplineNet netOf(const occ::handle<Geom_BSplineSurface> &);
gp_XY boundaryUV(const char *, const TopoDS_Face &, const TopoDS_Edge &, bool end);
int boundarySide(const char *, const TopoDS_Face &, const TopoDS_Edge &);
Constraint constraint(const char *, const TopoDS_Face &, const TopoDS_Edge &, Boundary, int);
// Derivative constraints only. Values are equation tolerances, scaled by the
// current metric; the defaults preserve the existing surface-tool contract.
struct JetFitSettings {
    double normalWeight = 1e-3, normalWithin = 1e-7;
    double curvatureWeight = 1e-6, curvatureWithin = 1e-9, metricFloor = 1;
};
void appendJetConditions(std::vector<FitCondition>&, const gp_XY&, const SurfaceJet& target,
                         const SurfaceJet& current, int order, const JetFitSettings& = {});
std::vector<FitCondition> conditions(const char *, const std::vector<Constraint> &,
                                     const occ::handle<Geom_BSplineSurface> &, double tolerance, bool normalFlow = false);
TopoDS_Face makeFace(const char *, const occ::handle<Geom_BSplineSurface> &, double tolerance,
                     TopAbs_Orientation orientation = TopAbs_FORWARD);
void measure(const char *, const TopoDS_Face &, const std::vector<Constraint> &, double tolerance, bool normalFlow = false);
} // namespace surfacing
EvolutionData alignCollinear(OcctKernel &, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, double, double, double,
                             double, double);
EvolutionData alignSurface(OcctKernel &, uint32_t, uint32_t, uint32_t, uint32_t, int, double, int = 0, int = 0, int = 0,
                           int = 0, int = 0, bool = false, double = 0, double = 1, bool = false, double = 0,
                           double = 0);
EvolutionData bridgeSurface(OcctKernel &, uint32_t, uint32_t, uint32_t, uint32_t, int, int, double, double,
                            double, int = 0);
EvolutionData squareSurface(OcctKernel &, const std::vector<uint32_t> &, const std::vector<uint32_t> &,
                            const std::vector<int> &, int, int, int, int, double, double = -1, bool = false);
std::vector<double> surfacePoles(OcctKernel &, uint32_t);
EvolutionData moveSurfacePoles(OcctKernel &, uint32_t, const std::vector<int> &, const std::vector<double> &,
                               const std::vector<uint32_t> &, const std::vector<uint32_t> &,
                               const std::vector<uint32_t> &, const std::vector<int> &, double);
EvolutionData raiseSurfaceDegree(OcctKernel &, uint32_t, int, int, double);
/// The B-spline that `faceSpline` data describes; fails naming `op` on malformed data.
occ::handle<Geom_BSplineSurface> splineFromData(const char *op, const std::vector<double> &);
} // namespace BRepExt
