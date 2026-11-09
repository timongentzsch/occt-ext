// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_CubicThrough.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>
#include <gp_Pln.hxx>
#include <gp_Ax3.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomConvert.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BSplCLib.hxx>
#include <Bnd_Box.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <array>
#include <functional>
#include <tuple>

namespace BRepExt
{
namespace
{
//! An edge's curve as a spline from 0 to 1, run the other way where `reversed`.
Handle(Geom_BSplineCurve) unitCurve(const TopoDS_Edge& edge, bool reversed)
{
  double                    a, b;
  const auto                curve = BRep_Tool::Curve(edge, a, b);
  Handle(Geom_BSplineCurve) spline =
    GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, a, b));
  if ((edge.Orientation() == TopAbs_REVERSED) != reversed)
    spline->Reverse();
  NCollection_Array1<double> knots(1, spline->NbKnots());
  spline->Knots(knots);
  BSplCLib::Reparametrize(0, 1, knots);
  spline->SetKnots(knots);
  return spline;
}

//! A guided skin between two profiles with G2 asked at an end: the sum of the surface straight
//! across between the guides, and the quintic in the guides' direction that has each profile, the
//! support's plane and bend along it, less the same between the guides' ends. Along a guide the
//! last two cancel, so each guide is kept; along a profile the first and last cancel. The two
//! supports give the plane and the bend, the guides how fast the surface leaves. Between flat
//! supports this is the guide swept straight across.
TopoDS_Face guidedSum(OcctKernel&                  kernel,
                      const char*                  op,
                      const std::vector<uint32_t>& edges,
                      uint32_t                     firstFace,
                      uint32_t                     lastFace,
                      int                          firstOrder,
                      int                          lastOrder,
                      double                       tolerance)
{
  const auto edgeAt  = [&](size_t i) { return TopoDS::Edge(get(kernel, edges[i])); };
  const auto startOf = [](const TopoDS_Edge& e) {
    double     a, b;
    const auto c = BRep_Tool::Curve(e, a, b);
    return c->Value(e.Orientation() == TopAbs_REVERSED ? b : a);
  };
  const auto endOf = [](const TopoDS_Edge& e) {
    double     a, b;
    const auto c = BRep_Tool::Curve(e, a, b);
    return c->Value(e.Orientation() == TopAbs_REVERSED ? a : b);
  };
  const auto first = edgeAt(0), last = edgeAt(1), g0 = edgeAt(2), g1 = edgeAt(3);
  // Each guide from the first profile to the last, and each profile from the first guide to the
  // second.
  const auto distanceTo = [](const gp_Pnt& p, const TopoDS_Edge& e) {
    BRepExtrema_DistShapeShape near(BRepBuilderAPI_MakeVertex(p).Vertex(), e);
    return near.Value();
  };
  const bool g0Back = distanceTo(startOf(g0), first) > distanceTo(endOf(g0), first);
  const bool g1Back = distanceTo(startOf(g1), first) > distanceTo(endOf(g1), first);
  const auto guideA = unitCurve(g0, g0Back), guideB = unitCurve(g1, g1Back);
  const bool firstBack =
    startOf(first).Distance(guideA->StartPoint()) > endOf(first).Distance(guideA->StartPoint());
  const bool lastBack =
    startOf(last).Distance(guideA->EndPoint()) > endOf(last).Distance(guideA->EndPoint());
  const auto rimA = unitCurve(first, firstBack), rimB = unitCurve(last, lastBack);
  // A profile at a share of its edge's parameter, which is where a support's jet is read; the
  // spline of a conic runs at another pace.
  const auto pointOn = [](const TopoDS_Edge& edge, bool back) {
    double                   a, b;
    const Handle(Geom_Curve) curve   = BRep_Tool::Curve(edge, a, b);
    const bool               flipped = (edge.Orientation() == TopAbs_REVERSED) != back;
    return [curve, a, b, flipped](double u) {
      return curve->Value(a + (b - a) * (flipped ? 1 - u : u));
    };
  };
  const auto rimPointA = pointOn(first, firstBack), rimPointB = pointOn(last, lastBack);
  const auto jetsOf = [&](uint32_t edgeId,
                          uint32_t faceId,
                          int      order,
                          bool     back) -> std::function<SurfaceJet(double)> {
    if (order < 1)
      return {};
    const auto boundary = surfacing::boundaryOf(kernel, op, edgeId, faceId);
    const bool flipped  = (boundary.edge.Orientation() == TopAbs_REVERSED) != back;
    return [boundary, flipped, op](double u) { return boundary.jet(op, flipped ? 1 - u : u); };
  };
  // How the surface leaves a profile at u: the guides' way there, in the support's plane, with
  // the support's bend in that direction.
  const auto leaving =
    [&](const std::function<SurfaceJet(double)>& jetAt, int order, double v, double u) {
      gp_Pnt ignored;
      gp_Vec a1, a2, b1, b2;
      guideA->D2(v, ignored, a1, a2);
      guideB->D2(v, ignored, b1, b2);
      gp_Vec d1 = a1 * (1 - u) + b1 * u, d2 = a2 * (1 - u) + b2 * u;
      if (order >= 1)
      {
        const auto jet = jetAt(u);
        d1 -= jet.n * d1.Dot(jet.n);
        if (order >= 2)
          d2 += jet.n * (jet.curvature(d1, d1) - d2.Dot(jet.n));
      }
      return std::pair{d1, d2};
    };
  const auto jetsA    = jetsOf(edges[0], firstFace, firstOrder, firstBack),
             jetsB    = jetsOf(edges[1], lastFace, lastOrder, lastBack);
  const auto splineOf = [](const std::vector<gp_Pnt>& values, const std::vector<double>& places) {
    const GeomAlgoExt::CubicThrough law(places, GeomAlgoExt::CubicThrough::Ends::Follow);
    const auto                      poles = law.through(values);
    NCollection_Array1<gp_Pnt>      array(1, law.poles);
    for (int k = 1; k <= law.poles; ++k)
      array(k) = poles[size_t(k - 1)];
    return Handle(Geom_BSplineCurve)(new Geom_BSplineCurve(array, law.knots, law.mults, 3));
  };
  // The fields along each profile as splines, and a rational profile the same way, for the sum
  // of poles, to the tolerance: over a span of 1 of the guides' parameter these are lengths.
  Handle(Geom_BSplineCurve) a1, a2, b1, b2, onA = rimA, onB = rimB;
  for (int n = 16;; n *= 2)
  {
    if (n > 4096)
      fail(op, no_solution, "the supports vary too much along the profiles to follow");
    std::vector<double> places;
    std::vector<gp_Pnt> va1, va2, vb1, vb2, pa, pb;
    for (int k = 0; k <= n; ++k)
    {
      const double u = double(k) / n;
      places.push_back(u);
      const auto [d1, d2] = leaving(jetsA, firstOrder, 0, u);
      const auto [e1, e2] = leaving(jetsB, lastOrder, 1, u);
      va1.emplace_back(d1.XYZ());
      va2.emplace_back(d2.XYZ());
      vb1.emplace_back(e1.XYZ());
      vb2.emplace_back(e2.XYZ());
      pa.push_back(rimPointA(u));
      pb.push_back(rimPointB(u));
    }
    a1 = splineOf(va1, places);
    a2 = splineOf(va2, places);
    b1 = splineOf(vb1, places);
    b2 = splineOf(vb2, places);
    if (rimA->IsRational())
      onA = splineOf(pa, places);
    if (rimB->IsRational())
      onB = splineOf(pb, places);
    double worst = 0;
    for (int k = 0; k < n; ++k)
    {
      const double u      = (k + 0.5) / n;
      const auto [d1, d2] = leaving(jetsA, firstOrder, 0, u);
      const auto [e1, e2] = leaving(jetsB, lastOrder, 1, u);
      worst               = std::max({worst,
                                      a1->Value(u).Distance(gp_Pnt(d1.XYZ())),
                                      a2->Value(u).Distance(gp_Pnt(d2.XYZ())),
                                      b1->Value(u).Distance(gp_Pnt(e1.XYZ())),
                                      b2->Value(u).Distance(gp_Pnt(e2.XYZ())),
                                      onA->Value(u).Distance(rimPointA(u)),
                                      onB->Value(u).Distance(rimPointB(u))});
    }
    if (worst <= tolerance / 4)
      break;
  }
  // The quintic: six rows of poles along the profiles.
  std::vector<Handle(Geom_BSplineCurve)> along{Handle(Geom_BSplineCurve)::DownCast(onA->Copy()),
                                               a1,
                                               a2,
                                               b1,
                                               b2,
                                               Handle(Geom_BSplineCurve)::DownCast(onB->Copy())};
  GeomAlgoExt::unify(along);
  const int                  nu = along[0]->NbPoles();
  NCollection_Array2<gp_Pnt> quintic(1, nu, 1, 6);
  for (int a = 1; a <= nu; ++a)
  {
    const gp_XYZ p0 = along[0]->Pole(a).XYZ(), d1 = along[1]->Pole(a).XYZ(),
                 d2 = along[2]->Pole(a).XYZ(), e1 = along[3]->Pole(a).XYZ(),
                 e2 = along[4]->Pole(a).XYZ(), p1 = along[5]->Pole(a).XYZ();
    quintic(a, 1) = gp_Pnt(p0);
    quintic(a, 2) = gp_Pnt(p0 + d1 / 5);
    quintic(a, 3) = gp_Pnt(p0 + d1 * 2 / 5 + d2 / 20);
    quintic(a, 4) = gp_Pnt(p1 - e1 * 2 / 5 + e2 / 20);
    quintic(a, 5) = gp_Pnt(p1 - e1 / 5);
    quintic(a, 6) = gp_Pnt(p1);
  }
  NCollection_Array1<double> bezierKnots(1, 2);
  NCollection_Array1<int>    quinticMults(1, 2), lineMults(1, 2);
  bezierKnots(1)  = 0;
  bezierKnots(2)  = 1;
  quinticMults(1) = quinticMults(2) = 6;
  lineMults(1) = lineMults(2)              = 2;
  Handle(Geom_BSplineSurface) profilesSkin = new Geom_BSplineSurface(quintic,
                                                                     along[0]->Knots(),
                                                                     bezierKnots,
                                                                     along[0]->Multiplicities(),
                                                                     quinticMults,
                                                                     along[0]->Degree(),
                                                                     5);
  // The same quintic between the ends of the profiles, straight across.
  NCollection_Array2<gp_Pnt> corners(1, 2, 1, 6);
  for (int b = 1; b <= 6; ++b)
  {
    corners(1, b) = quintic(1, b);
    corners(2, b) = quintic(nu, b);
  }
  Handle(Geom_BSplineSurface) crossings =
    new Geom_BSplineSurface(corners, bezierKnots, bezierKnots, lineMults, quinticMults, 1, 5);
  // The surface straight across between the guides.
  std::vector<Handle(Geom_BSplineCurve)> sides{guideA, guideB};
  GeomAlgoExt::unify(sides);
  const int                  nv = sides[0]->NbPoles();
  NCollection_Array2<gp_Pnt> straight(1, 2, 1, nv);
  for (int b = 1; b <= nv; ++b)
  {
    straight(1, b) = sides[0]->Pole(b);
    straight(2, b) = sides[1]->Pole(b);
  }
  Handle(Geom_BSplineSurface) guidesSkin = new Geom_BSplineSurface(straight,
                                                                   bezierKnots,
                                                                   sides[0]->Knots(),
                                                                   lineMults,
                                                                   sides[0]->Multiplicities(),
                                                                   1,
                                                                   sides[0]->Degree());
  // On one set of knots, and summed.
  const std::vector<Handle(Geom_BSplineSurface)> skins{profilesSkin, crossings, guidesSkin};
  int                                            du = 1, dv = 1;
  for (const auto& skin : skins)
  {
    du = std::max(du, skin->UDegree());
    dv = std::max(dv, skin->VDegree());
  }
  for (const auto& skin : skins)
    skin->IncreaseDegree(du, dv);
  for (const auto& from : skins)
    for (const auto& to : skins)
    {
      for (int k = 2; k < from->NbUKnots(); ++k)
        to->InsertUKnot(from->UKnot(k), from->UMultiplicity(k), 1e-12, false);
      for (int k = 2; k < from->NbVKnots(); ++k)
        to->InsertVKnot(from->VKnot(k), from->VMultiplicity(k), 1e-12, false);
    }
  for (int a = 1; a <= profilesSkin->NbUPoles(); ++a)
    for (int b = 1; b <= profilesSkin->NbVPoles(); ++b)
      profilesSkin->SetPole(a,
                            b,
                            gp_Pnt(profilesSkin->Pole(a, b).XYZ() + guidesSkin->Pole(a, b).XYZ()
                                   - crossings->Pole(a, b).XYZ()));
  const TopoDS_Face face = surfacing::makeFace(op, profilesSkin, tolerance);
  // What comes back is held to the supports apart from how it was made.
  for (const auto& [edgeId, faceId, order, v] :
       std::array{std::tuple{edges[0], firstFace, firstOrder, 0.0},
                  std::tuple{edges[1], lastFace, lastOrder, 1.0}})
  {
    if (order < 1)
      continue;
    const auto  reference = surfacing::boundaryOf(kernel, op, edgeId, faceId);
    TopoDS_Edge own;
    double      nearest = Precision::Infinite();
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
      const double      d = curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2)
                         .Distance(reference.point(0.5));
      if (d < nearest)
      {
        nearest = d;
        own     = TopoDS::Edge(it.Current());
      }
    }
    (void)v;
    surfacing::measure(op,
                       face,
                       {surfacing::constraint(op, face, own, reference, order)},
                       tolerance,
                       false);
  }
  return face;
}

//! A guided skin through three or more profiles with G1 or G2 asked at an end: the skin of the net
//! (skinGordon), plus a correction across each such end. At the end the correction turns the way
//! the skin leaves the profile into the support's tangent plane (a move along the support's normal
//! only) and, at G2, sets its bend across to the support's in that direction (again along the
//! normal), which is all that G1 and G2 ask of a surface that already has the profile as its edge.
//! Across, the correction is the quintic that has that way and bend at the end and is nothing, to
//! the second derivative, at the next profile, so every other profile is kept and the skin stays
//! C2 there. Along the profile it is the move less its straight interpolation between the guides,
//! so each guide is kept as well; where a guide leaves the end profile out of the support's plane
//! (or, at G2, with another bend), the surface follows the guide there, as between two profiles.
//! Null where the net's skin has weights other than one, to which a correction without weights
//! does not add; weights off one by rounding only are dropped.
TopoDS_Face guidedNet(OcctKernel&                  kernel,
                      const char*                  op,
                      const std::vector<uint32_t>& edges,
                      uint32_t                     firstFace,
                      uint32_t                     lastFace,
                      int                          firstOrder,
                      int                          lastOrder,
                      double                       tolerance)
{
  const size_t                last = edges.size() - 3;
  const std::vector<uint32_t> profiles(edges.begin(), edges.begin() + long(last) + 1),
    guides(edges.begin() + long(last) + 1, edges.end());
  const uint32_t netId = skinGordon(kernel, profiles, guides, false, tolerance).resultId;
  TopoDS_Face    net;
  for (TopExp_Explorer it(get(kernel, netId), TopAbs_FACE); it.More(); it.Next())
    net = TopoDS::Face(it.Current());
  release(kernel, netId);
  if (net.IsNull())
    return {};
  Handle(Geom_BSplineSurface) surface =
    Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(net)->Copy());
  if (surface.IsNull())
    return {};
  // Weights that differ from one by rounding only (GeomFill_Gordon's exact reparametrization
  // leaves some) move no pole by more than a hundredth of the tolerance when dropped.
  if (surface->IsURational() || surface->IsVRational())
  {
    double  spread = 0;
    Bnd_Box box;
    for (int a = 1; a <= surface->NbUPoles(); ++a)
      for (int b = 1; b <= surface->NbVPoles(); ++b)
      {
        spread = std::max(spread, std::abs(surface->Weight(a, b) - 1));
        box.Add(surface->Pole(a, b));
      }
    if (spread * std::sqrt(box.SquareExtent()) > tolerance / 100)
      return {};
    NCollection_Array2<gp_Pnt> poles(1, surface->NbUPoles(), 1, surface->NbVPoles());
    surface->Poles(poles);
    surface = new Geom_BSplineSurface(poles,
                                      surface->UKnots(),
                                      surface->VKnots(),
                                      surface->UMultiplicities(),
                                      surface->VMultiplicities(),
                                      surface->UDegree(),
                                      surface->VDegree());
  }
  // The first profile at v = 0 and the last at v = 1, both parameters from 0 to 1.
  const auto middleOf = [&](size_t i) {
    BRepAdaptor_Curve c(TopoDS::Edge(get(kernel, edges[i])));
    return c.Value((c.FirstParameter() + c.LastParameter()) / 2);
  };
  const auto sideNearest = [&](const gp_Pnt& p) {
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    const std::array<gp_Pnt, 4> sides{surface->Value(u0, (v0 + v1) / 2),
                                      surface->Value(u1, (v0 + v1) / 2),
                                      surface->Value((u0 + u1) / 2, v0),
                                      surface->Value((u0 + u1) / 2, v1)};
    int                         best = 0;
    for (int s = 1; s < 4; ++s)
      if (sides[size_t(s)].Distance(p) < sides[size_t(best)].Distance(p))
        best = s;
    return best;
  };
  int side = sideNearest(middleOf(0));
  if (side < 2)
  {
    surface->ExchangeUV();
    side += 2;
  }
  if (side == 3)
    surface->VReverse();
  if (sideNearest(middleOf(last)) != 3)
    return {};
  for (const bool alongU : {true, false})
  {
    NCollection_Array1<double> knots(1, alongU ? surface->NbUKnots() : surface->NbVKnots());
    alongU ? surface->UKnots(knots) : surface->VKnots(knots);
    BSplCLib::Reparametrize(0, 1, knots);
    alongU ? surface->SetUKnots(knots) : surface->SetVKnots(knots);
  }
  // How far across the correction may reach: to the nearest place of the next profile.
  const auto reach = [&](size_t i, bool fromStart) {
    BRepAdaptor_Curve     c(TopoDS::Edge(get(kernel, edges[i])));
    double                bound = fromStart ? 1.0 : 0.0;
    ShapeAnalysis_Surface on(surface);
    for (int s = 0; s <= 16; ++s)
    {
      const double v = on.ValueOfUV(c.Value(c.FirstParameter()
                                            + (c.LastParameter() - c.FirstParameter()) * s / 16),
                                    Precision::Confusion())
                         .Y();
      bound = fromStart ? std::min(bound, v) : std::max(bound, v);
    }
    return fromStart ? bound : 1 - bound;
  };
  std::vector<Handle(Geom_BSplineSurface)> corrections;
  for (const auto& [edgeId, faceId, order, atStart] :
       std::array{std::tuple{edges[0], firstFace, firstOrder, true},
                  std::tuple{edges[last], lastFace, lastOrder, false}})
  {
    if (order < 1)
      continue;
    const double width = reach(atStart ? 1 : last - 1, atStart);
    if (!(width > Precision::PConfusion()))
      return {};
    const auto   boundary = surfacing::boundaryOf(kernel, op, edgeId, faceId);
    const double v        = atStart ? 0 : 1;
    const int    asked    = order; // a lambda captures no structured binding before C++20
    // The move of the way out and of the bend at u, both along the support's normal there.
    const auto moves = [&](double u) {
      gp_Pnt p;
      gp_Vec su, sv, suu, svv, suv;
      surface->D2(u, v, p, su, sv, suu, svv, suv);
      GeomAPI_ProjectPointOnCurve onRim(p, boundary.curve, boundary.first, boundary.last);
      const double f   = onRim.NbPoints() > 0 ? onRim.LowerDistanceParameter() : boundary.first;
      const auto   jet = boundary.jet(op, (f - boundary.first) / (boundary.last - boundary.first));
      const gp_Vec d1  = sv - jet.n * sv.Dot(jet.n);
      const gp_Vec d2  = asked >= 2 ? jet.n * (jet.curvature(d1, d1) - svv.Dot(jet.n)) : gp_Vec();
      return std::pair{d1 - sv, d2};
    };
    const auto [m1a, m2a] = moves(0);
    const auto [m1b, m2b] = moves(1);
    const auto corrected  = [&](double u) {
      const auto [m1, m2] = moves(u);
      return std::pair{m1 - m1a * (1 - u) - m1b * u, m2 - m2a * (1 - u) - m2b * u};
    };
    const auto splineOf = [](const std::vector<gp_Pnt>& values, const std::vector<double>& places) {
      const GeomAlgoExt::CubicThrough law(places, GeomAlgoExt::CubicThrough::Ends::Follow);
      const auto                      poles = law.through(values);
      NCollection_Array1<gp_Pnt>      array(1, law.poles);
      for (int k = 1; k <= law.poles; ++k)
        array(k) = poles[size_t(k - 1)];
      return Handle(Geom_BSplineCurve)(new Geom_BSplineCurve(array, law.knots, law.mults, 3));
    };
    // The two moves as splines along the profile, to a quarter of the tolerance: on the unit
    // parameter across they are lengths.
    Handle(Geom_BSplineCurve) way, bend;
    for (int n = 16;; n *= 2)
    {
      if (n > 4096)
        fail(op, no_solution, "the supports vary too much along the profiles to follow");
      std::vector<double> places;
      std::vector<gp_Pnt> ways, bends;
      for (int k = 0; k <= n; ++k)
      {
        const auto [m1, m2] = corrected(double(k) / n);
        places.push_back(double(k) / n);
        ways.emplace_back(m1.XYZ());
        bends.emplace_back(m2.XYZ());
      }
      way          = splineOf(ways, places);
      bend         = splineOf(bends, places);
      double worst = 0;
      for (int k = 0; k < n; ++k)
      {
        const auto [m1, m2] = corrected((k + 0.5) / n);
        worst               = std::max({worst,
                                        way->Value((k + 0.5) / n).Distance(gp_Pnt(m1.XYZ())),
                                        bend->Value((k + 0.5) / n).Distance(gp_Pnt(m2.XYZ()))});
      }
      if (worst <= tolerance / 4)
        break;
    }
    // The quintic across, on [0, width] from the end: value nothing, first derivative `way`,
    // second `bend` at the end, and nothing to the second derivative at `width`; nothing beyond.
    std::vector<Handle(Geom_BSplineCurve)> along{way, bend};
    GeomAlgoExt::unify(along);
    const int                  nu = along[0]->NbPoles();
    NCollection_Array2<gp_Pnt> poles(1, nu, 1, 11);
    for (int a = 1; a <= nu; ++a)
    {
      const gp_XYZ d1 = along[0]->Pole(a).XYZ(), d2 = along[1]->Pole(a).XYZ();
      for (int b = 1; b <= 11; ++b)
        poles(a, b) = gp_Pnt(0, 0, 0);
      // The Bezier poles of the quintic on its stretch (two Bezier pieces sharing the sixth
      // pole): from the end, nothing, then the way out over a fifth of the stretch, then twice
      // that with the bend; the other three nothing.
      const double sign          = atStart ? 1 : -1;
      poles(a, atStart ? 2 : 10) = gp_Pnt(d1 * (sign * width / 5));
      poles(a, atStart ? 3 : 9)  = gp_Pnt(d1 * (sign * 2 * width / 5) + d2 * (width * width / 20));
    }
    NCollection_Array1<double> vKnots(1, 3);
    NCollection_Array1<int>    vMults(1, 3);
    vKnots(1) = 0;
    vKnots(2) = atStart ? width : 1 - width;
    vKnots(3) = 1;
    vMults(1) = vMults(3)                  = 6;
    vMults(2)                              = 5;
    Handle(Geom_BSplineSurface) correction = new Geom_BSplineSurface(poles,
                                                                     along[0]->Knots(),
                                                                     vKnots,
                                                                     along[0]->Multiplicities(),
                                                                     vMults,
                                                                     3,
                                                                     5);
    // Nothing to the second derivative on both sides of the stretch's end: C2 there.
    if (!correction->RemoveVKnot(2, 3, Precision::Confusion()))
      return {};
    corrections.push_back(correction);
  }
  // On one set of knots, and summed.
  std::vector<Handle(Geom_BSplineSurface)> all{surface};
  all.insert(all.end(), corrections.begin(), corrections.end());
  int du = 1, dv = 1;
  for (const auto& s : all)
  {
    du = std::max(du, s->UDegree());
    dv = std::max(dv, s->VDegree());
  }
  for (const auto& s : all)
    s->IncreaseDegree(du, dv);
  for (const auto& from : all)
    for (const auto& to : all)
    {
      for (int k = 2; k < from->NbUKnots(); ++k)
        to->InsertUKnot(from->UKnot(k), from->UMultiplicity(k), 1e-12, false);
      for (int k = 2; k < from->NbVKnots(); ++k)
        to->InsertVKnot(from->VKnot(k), from->VMultiplicity(k), 1e-12, false);
    }
  for (const auto& correction : corrections)
    for (int a = 1; a <= surface->NbUPoles(); ++a)
      for (int b = 1; b <= surface->NbVPoles(); ++b)
        surface->SetPole(a, b, gp_Pnt(surface->Pole(a, b).XYZ() + correction->Pole(a, b).XYZ()));
  const TopoDS_Face face = surfacing::makeFace(op, surface, tolerance);
  // What comes back is held to the supports apart from how it was made.
  for (const auto& [edgeId, faceId, order] :
       std::array{std::tuple{edges[0], firstFace, firstOrder},
                  std::tuple{edges[last], lastFace, lastOrder}})
  {
    if (order < 1)
      continue;
    const auto  reference = surfacing::boundaryOf(kernel, op, edgeId, faceId);
    TopoDS_Edge own;
    double      nearest = Precision::Infinite();
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
      const double      d = curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2)
                         .Distance(reference.point(0.5));
      if (d < nearest)
      {
        nearest = d;
        own     = TopoDS::Edge(it.Current());
      }
    }
    surfacing::measure(op,
                       face,
                       {surfacing::constraint(op, face, own, reference, order)},
                       tolerance,
                       false);
  }
  return face;
}
} // namespace

EvolutionData skinGuides(OcctKernel&                  kernel,
                         const std::vector<uint32_t>& edges,
                         uint32_t                     firstFace,
                         uint32_t                     lastFace,
                         int                          firstOrder,
                         int                          lastOrder,
                         double                       tolerance)
{
  const char* op = "skinGuides";
  if (edges.size() < 4 || firstOrder < 0 || firstOrder > 2 || lastOrder < 0 || lastOrder > 2
      || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op,
         invalid_input,
         "at least two profiles, two guides and G0-G2 end constraints are required");
  BRepOffsetAPI_MakeFilling
                          fill(3, 50, 5, false, tolerance / 10, tolerance / 10, 1e-5, 1e-4, 8, 128);
  BRepBuilderAPI_MakeWire boundary;
  // Connectivity is checked independently from filling, which otherwise accepts
  // disconnected constraints.
  for (auto id : edges)
  {
    if (get(kernel, id).ShapeType() != TopAbs_EDGE)
      fail(op, invalid_input, "each profile and guide must be a single edge");
  }
  const int last = int(edges.size()) - 3;
  for (int i : {0, last + 1, last, last + 2})
    boundary.Add(TopoDS::Edge(get(kernel, edges[i])));
  if (!boundary.IsDone() || !boundary.Wire().Closed())
    fail(op, invalid_input, "profiles and guides must form one closed four-sided boundary");
  if (edges.size() == 4 && (firstOrder == 2 || lastOrder == 2))
  {
    if ((firstOrder > 0 && !firstFace) || (lastOrder > 0 && !lastFace))
      fail(op, invalid_input, "G1/G2 requires an adjacent support face");
    return stored(
      kernel,
      guidedSum(kernel, op, edges, firstFace, lastFace, firstOrder, lastOrder, tolerance));
  }
  for (int i = 0; i < int(edges.size()); ++i)
  {
    const auto     edge    = TopoDS::Edge(get(kernel, edges[i]));
    const int      order   = i == 0 ? firstOrder : i == last ? lastOrder : 0;
    const uint32_t support = i == 0 ? firstFace : i == last ? lastFace : 0;
    if (order)
    {
      if (!support || get(kernel, support).ShapeType() != TopAbs_FACE)
        fail(op, invalid_input, "G1/G2 requires an adjacent support face");
      subShapes(kernel, op, get(kernel, support), {edges[i]}, TopAbs_EDGE);
      // This OCCT build forwards GeomAbs_Shape numerically to the plate
      // constraint, which accepts integer orders 0, 1, 2 (not enum G2=3).
      fill.Add(edge, TopoDS::Face(get(kernel, support)), static_cast<GeomAbs_Shape>(order));
    }
    else
      fill.Add(edge, GeomAbs_C0, i == 0 || i >= last);
  }
  // Coplanar inputs admit the zero-bending plate solution. Keep its existing
  // carrier and trims: a network's unit rectangle is a different representation
  // of the same plane. This construction also retains every interior profile.
  // The plane search uses the requested model tolerance; the plate retains all
  // its original physical acceptance checks below.
  TopoDS_Compound planeInputs;
  BRep_Builder planeBuilder;
  planeBuilder.MakeCompound(planeInputs);
  for (const auto id : edges)
    planeBuilder.Add(planeInputs, get(kernel, id));
  bool planarSupports = true;
  for (const auto id : {firstOrder ? firstFace : 0u, lastOrder ? lastFace : 0u})
  {
    if (!id)
      continue;
    const auto face = TopoDS::Face(get(kernel, id));
    planarSupports = planarSupports && BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane;
    planeBuilder.Add(planeInputs, face);
  }
  const bool planar = planarSupports && BRepBuilderAPI_FindPlane(planeInputs, tolerance).Found();
  bool attemptedNetwork = false;
  const auto network = [&]() -> TopoDS_Face {
    attemptedNetwork = true;
    try
    {
      return guidedNet(kernel, op, edges, firstFace, lastFace, firstOrder, lastOrder, tolerance);
    }
    catch (const Failure&)
    {
    }
    catch (const Standard_Failure&)
    {
    }
    return {};
  };
  // Curved compatible networks avoid the dense plate solve. Coplanar networks
  // retain the plate and its existing recovery; each construction runs at most once.
  if (edges.size() > 4 && !planar)
  {
    const TopoDS_Face face = network();
    if (!face.IsNull())
      return stored(kernel, face);
  }
  try
  {
    fill.Build();
    if (!fill.IsDone() || fill.G0Error() > tolerance || !BRepCheck_Analyzer(fill.Shape()).IsValid())
      fail(op, no_solution, "guided skin did not meet its boundary tolerance");
    if ((firstOrder > 0 || lastOrder > 0) && fill.G1Error() > 1e-4)
      fail(op, no_solution, "guided skin did not meet G1");
    if ((firstOrder > 1 || lastOrder > 1) && fill.G2Error() > 1e-3)
      fail(op, no_solution, "guided skin did not meet G2: " + std::to_string(fill.G2Error()));
    return stored(kernel, fill.Shape());
  }
  catch (const Failure&)
  {
    if (edges.size() > 4 && !attemptedNetwork && (firstOrder > 0 || lastOrder > 0))
    {
      const TopoDS_Face face = network();
      if (!face.IsNull())
        return stored(kernel, face);
    }
    throw;
  }
}
} // namespace BRepExt

#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array2.hxx>

namespace BRepExt
{
EvolutionData skinPeriodic(OcctKernel& kernel, const std::vector<uint32_t>& ids, double tolerance)
{
  const char* op = "skinPeriodic";
  if (ids.size() < 3 || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "at least three sections and positive tolerance required");
  constexpr int              samples = 65;
  NCollection_Array2<gp_Pnt> grid(1, int(ids.size()), 1, samples);
  std::vector<TopoDS_Edge>   edges;
  for (size_t i = 0; i < ids.size(); ++i)
  {
    const auto&     shape = get(kernel, ids[i]);
    TopExp_Explorer it(shape, TopAbs_EDGE);
    if (!it.More())
      fail(op, invalid_input, "empty profile");
    auto edge = TopoDS::Edge(it.Current());
    it.Next();
    if (it.More() || edge.Closed())
      fail(op, unsupported_geometry, "periodic skin requires single open-edge sections");
    BRepAdaptor_Curve curve(edge);
    double            a = curve.FirstParameter(), b = curve.LastParameter();
    if (edge.Orientation() == TopAbs_REVERSED)
      std::swap(a, b);
    if (i > 0)
    {
      const auto   previousStart = grid(int(i), 1), previousEnd = grid(int(i), samples);
      const double same =
        previousStart.SquareDistance(curve.Value(a)) + previousEnd.SquareDistance(curve.Value(b));
      const double reverse =
        previousStart.SquareDistance(curve.Value(b)) + previousEnd.SquareDistance(curve.Value(a));
      if (reverse < same)
      {
        std::swap(a, b);
        edge.Reverse();
      }
    }
    edges.push_back(edge);
    for (int j = 0; j < samples; ++j)
      grid(int(i) + 1, j + 1) = curve.Value(a + (b - a) * j / (samples - 1));
  }
  GeomAPI_PointsToBSplineSurface fit;
  fit.Interpolate(grid, Approx_IsoParametric, true);
  if (!fit.IsDone())
    fail(op, no_solution, "periodic interpolation failed");
  auto surface = fit.Surface();
  if (!surface->IsUPeriodic())
    fail(op, no_solution, "interpolation did not close tangentially");
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  // Check between interpolation sites as well. This is a sampled tolerance
  // check, not a certified error bound for arbitrarily oscillatory profiles.
  for (size_t i = 0; i < edges.size(); ++i)
  {
    BRepAdaptor_Curve curve(edges[i]);
    double            a = curve.FirstParameter(), b = curve.LastParameter();
    if (edges[i].Orientation() == TopAbs_REVERSED)
      std::swap(a, b);
    for (int j = 0; j < 2 * samples - 1; ++j)
    {
      const double f        = double(j) / (2 * samples - 2);
      const auto   expected = curve.Value(a + (b - a) * f);
      if (expected.Distance(surface->Value(u0 + (u1 - u0) * i / ids.size(), v0 + (v1 - v0) * f))
          > tolerance)
        fail(op,
             no_solution,
             "periodic skin misses a section between samples; rebuild or simplify the section");
    }
  }
  auto face = BRepBuilderAPI_MakeFace(surface, tolerance).Face();
  if (!BRepCheck_Analyzer(face).IsValid())
    fail(op, no_solution, "invalid periodic surface");
  return stored(kernel, face);
}
} // namespace BRepExt
