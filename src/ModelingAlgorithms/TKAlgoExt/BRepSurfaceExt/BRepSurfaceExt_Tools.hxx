// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

namespace BRepExt
{
namespace surfacing
{
struct Boundary
{
  TopoDS_Edge             edge;
  TopoDS_Face             face;
  occ::handle<Geom_Curve> curve;
  double                  first = 0, last = 0;
  gp_Pnt                  point(double f) const;
  SurfaceJet              jet(const char* op, double f) const;
};

struct Constraint
{
  Boundary reference;
  gp_XY    start, end;
  int      order = 0;
  // How far the surface's way off the boundary leans along it at its two ends, as the tangent of
  // the angle from straight across: towards the other end when positive. Between them it runs
  // evenly. Only where `skewed`. Unsourced convention: the origin and sense of the angle are
  // chosen, not derived.
  bool   skewed    = false;
  double skewStart = 0, skewEnd = 0;

  double skewAt(double f) const { return skewStart * (1 - f) - skewEnd * f; }
};

//! The face behind `face`. Fails naming `op` where it is not a face.
Standard_EXPORT TopoDS_Face faceOf(OcctKernel& kernel, const char* op, uint32_t face);
//! `edge` of `face` with its curve and range. Fails where the edge is not on the face.
Standard_EXPORT Boundary boundaryOf(OcctKernel& kernel,
                                    const char* op,
                                    uint32_t    edge,
                                    uint32_t    face);
//! A copy of the surface of an untrimmed rectangular, non-periodic B-spline face; with `fitting`,
//! one that is not rational and of degree 5 at most. Fails naming `op` otherwise.
Standard_EXPORT occ::handle<Geom_BSplineSurface> editable(const char*        op,
                                                          const TopoDS_Face& face,
                                                          bool               fitting = true);
//! The degrees, knot sequences and poles of `surface`, U-major.
Standard_EXPORT SplineNet netOf(const occ::handle<Geom_BSplineSurface>& surface);
//! The parameters of the start, or with `end` the end, of `edge` on `face`.
Standard_EXPORT gp_XY boundaryUV(const char*        op,
                                 const TopoDS_Face& face,
                                 const TopoDS_Edge& edge,
                                 bool               end);
//! The side of the parameter rectangle that `edge` lies on whole: 0 u-min, 1 u-max, 2 v-min,
//! 3 v-max. Fails naming `op` otherwise.
Standard_EXPORT int boundarySide(const char* op, const TopoDS_Face& face, const TopoDS_Edge& edge);
//! The fit target that makes `edge` of `face` meet `reference` to `order` 0, 1 or 2 (G0 to G2).
Standard_EXPORT Constraint constraint(const char*        op,
                                      const TopoDS_Face& face,
                                      const TopoDS_Edge& edge,
                                      Boundary           reference,
                                      int                order);

//! The fit conditions (positions, normals, curvatures, skew) sampled along each constraint.
Standard_EXPORT std::vector<FitCondition> conditions(
  const char*                             op,
  const std::vector<Constraint>&          constraints,
  const occ::handle<Geom_BSplineSurface>& surface,
  double                                  tolerance,
  bool                                    normalFlow = false);
//! A valid face on `surface`. Fails naming `op` where it is not valid.
Standard_EXPORT TopoDS_Face makeFace(const char*                             op,
                                     const occ::handle<Geom_BSplineSurface>& surface,
                                     double                                  tolerance,
                                     TopAbs_Orientation orientation          = TopAbs_FORWARD,
                                     bool               useTopologyTolerance = false);
//! Check a fitted face against its constraints. Fails naming the gap, angle or curvature missed.
Standard_EXPORT void measure(const char*                    op,
                             const TopoDS_Face&             face,
                             const std::vector<Constraint>& constraints,
                             double                         tolerance,
                             bool                           normalFlow = false);
} // namespace surfacing

//! The shape of a bridge between two boundary edges.
enum class BridgeShape
{
  Smooth  = 0, //!< to the orders asked at each edge
  Chamfer = 1  //!< G0 at both edges
};

//! The row of poles next to `edge` of `face` put in line with the hull of `reference` on
//! `support`; with `symmetric`, both faces change. Skews in radians. Fails where the rows fold.
Standard_EXPORT EvolutionData alignCollinear(OcctKernel& kernel,
                                             uint32_t    face,
                                             uint32_t    edge,
                                             uint32_t    reference,
                                             uint32_t    support,
                                             bool        symmetric,
                                             bool        matchRatio,
                                             double      partialStart,
                                             double      partialEnd,
                                             double      skewStart,
                                             double      skewEnd,
                                             double      tolerance);
//! `edge` of a spline `face` made to meet `reference` on `support` to `order` 0 to 2 (G0 to G2),
//! raising degrees (to 5) and spans (to 9) exactly. Skews in radians. Fails where either would
//! fall.
Standard_EXPORT EvolutionData alignSurface(OcctKernel& kernel,
                                           uint32_t    face,
                                           uint32_t    edge,
                                           uint32_t    reference,
                                           uint32_t    support,
                                           int         order,
                                           double      tolerance,
                                           int         extraRows           = 0,
                                           int         uDegree             = 0,
                                           int         vDegree             = 0,
                                           int         uSpans              = 0,
                                           int         vSpans              = 0,
                                           bool        normalFlow          = false,
                                           double      partialStart        = 0,
                                           double      partialEnd          = 1,
                                           bool        skewed              = false,
                                           double      skewStart           = 0,
                                           double      skewEnd             = 0,
                                           double      inputShapeInfluence = 0);
//! A face bridging two boundary edges, to the orders 0 to 2 at each with tensions in (0, 10].
//! Fails where the edges are one.
Standard_EXPORT EvolutionData bridgeSurface(OcctKernel& kernel,
                                            uint32_t    firstEdge,
                                            uint32_t    firstFace,
                                            uint32_t    secondEdge,
                                            uint32_t    secondFace,
                                            int         orderA,
                                            int         orderB,
                                            double      tensionA,
                                            double      tensionB,
                                            double      tolerance,
                                            BridgeShape shape = BridgeShape::Smooth);
//! A four-sided spline face from two to four boundary `edges`, the missing sides made, each
//! smooth to its support to its order 0 to 2. Fails where they form no quad.
Standard_EXPORT EvolutionData squareSurface(OcctKernel&                  kernel,
                                            const std::vector<uint32_t>& edges,
                                            const std::vector<uint32_t>& supports,
                                            const std::vector<int>&      orders,
                                            int                          uDegree,
                                            int                          vDegree,
                                            int                          uSpans,
                                            int                          vSpans,
                                            double                       tolerance,
                                            double                       flatness   = -1,
                                            bool                         normalFlow = false);
//! The pole counts in u and v, each pole as x, y, z, weight (U-major), then the unit normal at
//! each pole's Greville point. Fails on a face that is not a B-spline.
Standard_EXPORT std::vector<double> surfacePoles(OcctKernel& kernel, uint32_t face);
//! The poles at the (u, v) `indices` from 0 moved by `offsets` (x, y, z each), the joins to up to
//! four neighbours kept to their `orders`. Fails above 400 poles.
Standard_EXPORT EvolutionData moveSurfacePoles(OcctKernel&                  kernel,
                                               uint32_t                     face,
                                               const std::vector<int>&      indices,
                                               const std::vector<double>&   offsets,
                                               const std::vector<uint32_t>& edges,
                                               const std::vector<uint32_t>& references,
                                               const std::vector<uint32_t>& supports,
                                               const std::vector<int>&      orders,
                                               double                       tolerance);
//! A B-spline face with its degrees raised by `deltaU` and `deltaV`. Fails past degree 25.
Standard_EXPORT EvolutionData
  raiseSurfaceDegree(OcctKernel& kernel, uint32_t face, int deltaU, int deltaV, double tolerance);
//! The B-spline that `faceSpline` data describes; fails naming `op` on malformed data.
Standard_EXPORT occ::handle<Geom_BSplineSurface> splineFromData(const char*                op,
                                                                const std::vector<double>& data);
} // namespace BRepExt
