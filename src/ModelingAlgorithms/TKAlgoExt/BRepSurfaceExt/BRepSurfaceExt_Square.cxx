// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_UniformKnots.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt
{
namespace
{
double greville(const std::vector<double>& k, int p, int i)
{
  double sum = 0;
  for (int j = 1; j <= p; ++j)
    sum += k[i + j];
  return sum / p;
}
} // namespace

EvolutionData squareSurface(OcctKernel&                  k,
                            const std::vector<uint32_t>& inputEdges,
                            const std::vector<uint32_t>& inputSupports,
                            const std::vector<int>&      inputOrders,
                            int                          pu,
                            int                          pv,
                            int                          su,
                            int                          sv,
                            double                       tolerance,
                            double                       flatness,
                            bool                         normalFlow)
{
  const char* op = "squareSurface";
  if (inputEdges.size() < 2 || inputEdges.size() > 4 || inputSupports.size() != inputEdges.size()
      || inputOrders.size() != inputEdges.size() || pu < 1 || pv < 1 || pu > 5 || pv > 5 || su < 1
      || sv < 1 || su > 9 || sv > 9 || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op,
         invalid_input,
         "two to four boundaries, degrees 1..5, spans 1..9 and a positive finite tolerance are "
         "required");
  if (!std::isfinite(flatness) || (flatness != -1 && (flatness < 0 || flatness > 1)))
    fail(op, invalid_input, "Flatness must be in 0..1");
  auto                  edges = inputEdges, supports = inputSupports;
  auto                  orders = inputOrders;
  std::vector<uint32_t> generated;

  // Generated boundaries are local arena entries, released on every exit.
  struct Release
  {
    OcctKernel&            k;
    std::vector<uint32_t>& ids;

    ~Release()
    {
      for (auto id : ids)
        BRepExt::release(k, id);
    }
  } release{k, generated};

  const auto curve = [&](uint32_t id) {
    if (get(k, id).ShapeType() != TopAbs_EDGE)
      fail(op, invalid_input, "boundaries must be edges");
    double a, b;
    auto   c = BRep_Tool::Curve(TopoDS::Edge(get(k, id)), a, b);
    if (c.IsNull())
      fail(op, invalid_input, "boundaries need spatial curves");
    return occ::handle<Geom_Curve>(new Geom_TrimmedCurve(c, a, b));
  };
  const auto add = [&](const TopoDS_Edge& edge) {
    const auto id = store(k, edge);
    generated.push_back(id);
    edges.push_back(id);
    supports.push_back(0);
    orders.push_back(0);
  };
  if (edges.size() == 2)
  {
    const auto   a = curve(edges[0]), b = curve(edges[1]);
    const gp_Pnt endsA[] = {a->Value(a->FirstParameter()), a->Value(a->LastParameter())},
                 endsB[] = {b->Value(b->FirstParameter()), b->Value(b->LastParameter())};
    int nearA = 0, nearB = 0, touches = 0;
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        if (endsA[i].Distance(endsB[j]) <= tolerance)
        {
          nearA = i;
          nearB = j;
          ++touches;
        }
    if (touches > 1)
      fail(op, unsupported_geometry, "two boundaries sharing both ends do not define a quad");
    if (touches)
    {
      // A translation sweep: the two selected sides meet, and copies of each
      // side through the other's far endpoint close the opposite sides.
      add(BRepBuilderAPI_MakeEdge(
            occ::down_cast<Geom_Curve>(a->Translated(gp_Vec(endsB[nearB], endsB[1 - nearB]))))
            .Edge());
      add(BRepBuilderAPI_MakeEdge(
            occ::down_cast<Geom_Curve>(b->Translated(gp_Vec(endsA[nearA], endsA[1 - nearA]))))
            .Edge());
    }
    else
    {
      const bool reverse = endsA[0].Distance(endsB[1]) + endsA[1].Distance(endsB[0])
                           < endsA[0].Distance(endsB[0]) + endsA[1].Distance(endsB[1]);
      add(BRepBuilderAPI_MakeEdge(endsA[0], endsB[reverse ? 1 : 0]).Edge());
      add(BRepBuilderAPI_MakeEdge(endsA[1], endsB[reverse ? 0 : 1]).Edge());
    }
  }
  else if (edges.size() == 3)
  {
    // Three touching sides in a chain. The missing fourth side is straight;
    // a closed triangle has no fourth corner and is refused rather than collapsed.
    std::vector<gp_Pnt> ends;
    for (auto id : edges)
    {
      const auto c = curve(id);
      ends.push_back(c->Value(c->FirstParameter()));
      ends.push_back(c->Value(c->LastParameter()));
    }
    std::vector<gp_Pnt> loose;
    for (size_t i = 0; i < ends.size(); ++i)
    {
      int neighbours = 0;
      for (size_t j = 0; j < ends.size(); ++j)
        if (i != j && ends[i].Distance(ends[j]) <= tolerance)
          ++neighbours;
      if (!neighbours)
        loose.push_back(ends[i]);
      if (neighbours > 1)
        fail(op, unsupported_geometry, "three boundaries must form an unbranched open chain");
    }
    if (loose.size() != 2)
      fail(op, unsupported_geometry, "three boundaries must form an open chain");
    add(BRepBuilderAPI_MakeEdge(loose[0], loose[1]).Edge());
  }
  // Ask Patch for its four-sided Coons path, without forcing these explicit knots on it.
  const auto seed     = patch(k,
                          edges,
                          std::vector<uint32_t>(4, 0),
                          std::vector<int>(4, 0),
                              {},
                              {},
                          tolerance,
                          0.1,
                          0.05,
                          1);
  const auto seedFace = surfacing::faceOf(k, op, seed.evolution.resultId);
  BRepExt::release(k, seed.evolution.resultId);
  if (seed.values[0] != 1)
    fail(op,
         unsupported_geometry,
         "the loop does not admit Patch's untrimmed four-sided construction");
  const auto original = surfacing::editable(op, seedFace, false);
  SplineNet  net{pu, pv, GeomAlgoExt::uniformKnots(pu, su), GeomAlgoExt::uniformKnots(pv, sv), {}};
  double     a, b, c, d;
  original->Bounds(a, b, c, d);
  for (int i = 0; i < net.countU(); ++i)
    for (int j = 0; j < net.countV(); ++j)
      net.poles.push_back(original
                            ->Value(a + (b - a) * greville(net.knotsU, pu, i),
                                    c + (d - c) * greville(net.knotsV, pv, j))
                            .XYZ());
  auto                               surface = surfaceOf(net);
  std::vector<surfacing::Constraint> targets;
  // Match each original boundary to its side of the Coons surface, irrespective of input order.
  for (size_t i = 0; i < 4; ++i)
  {
    if (orders[i] < -1 || orders[i] > 2 || (orders[i] > 0 && !supports[i]))
      fail(op, invalid_input, "G1/G2 need a support face");
    const auto& shape = get(k, edges[i]);
    if (shape.ShapeType() != TopAbs_EDGE)
      fail(op, invalid_input, "boundaries must be edges");
    auto                edge = TopoDS::Edge(shape);
    double              first, last;
    auto                curve = BRep_Tool::Curve(edge, first, last);
    surfacing::Boundary ref{edge, {}, curve, first, last};
    if (supports[i])
      ref = surfacing::boundaryOf(k, op, edges[i], supports[i]);
    const auto  p       = ref.point(.5);
    int         side    = 0;
    double      nearest = HUGE_VAL;
    const gp_XY mids[]  = {{0, .5}, {1, .5}, {.5, 0}, {.5, 1}};
    for (int s = 0; s < 4; ++s)
    {
      const auto q = original->Value(a + (b - a) * mids[s].X(), c + (d - c) * mids[s].Y());
      if (q.Distance(p) < nearest)
      {
        nearest = q.Distance(p);
        side    = s;
      }
    }
    gp_XY start   = side < 2 ? gp_XY(side, 0) : gp_XY(0, side - 2),
          end     = side < 2 ? gp_XY(side, 1) : gp_XY(1, side - 2);
    const auto p0 = original->Value(a + (b - a) * start.X(), c + (d - c) * start.Y());
    const auto p1 = original->Value(a + (b - a) * end.X(), c + (d - c) * end.Y());
    if (p0.Distance(ref.point(1)) + p1.Distance(ref.point(0))
        < p0.Distance(ref.point(0)) + p1.Distance(ref.point(1)))
      std::swap(ref.first, ref.last);
    targets.push_back({ref, start, end, orders[i]});
  }
  std::vector<char> fixed(net.poles.size(), 0);
  // Refit the Coons seed as well as its rim. Boundary jets alone leave the interior
  // underdetermined; these soft wishes retain its shape while the hard jets win.
  std::vector<FitCondition> seedConditions;
  for (int i = 0; i <= 16; ++i)
    for (int j = 0; j <= 16; ++j)
    {
      const double u = double(i) / 16, v = double(j) / 16;
      const auto   p = original->Value(a + (b - a) * u, c + (d - c) * v);
      for (int axis = 1; axis <= 3; ++axis)
      {
        gp_XYZ direction(0, 0, 0);
        direction.SetCoord(axis, 1);
        seedConditions.push_back({u, v, {{0, 0, direction}}, 1, p.Coord(axis), 1e-5, true});
      }
    }
  std::vector<surfacing::Constraint> held;
  std::copy_if(targets.begin(), targets.end(), std::back_inserter(held), [](const auto& target) {
    return target.order >= 0;
  });
  const FitEnergy energy{(flatness < 0 ? 1 : 1 - flatness) * 1e-11,
                         (flatness < 0 ? .01 : flatness) * 1e-11};
  const bool g2 = std::any_of(orders.begin(), orders.end(), [](int order) { return order == 2; });
  for (int pass = 0; pass < (g2 ? 8 : 1); ++pass)
  {
    auto constraints = surfacing::conditions(op, held, surface, tolerance, normalFlow);
    constraints.insert(constraints.end(), seedConditions.begin(), seedConditions.end());
    fitNet(op, net, fixed, constraints, energy, {});
    surface = surfaceOf(net);
  }
  auto       orientation = TopAbs_FORWARD;
  const auto oriented    = std::find_if(targets.begin(), targets.end(), [](const auto& target) {
    return target.order > 0;
  });
  if (oriented != targets.end())
  {
    const auto uv = (oriented->start + oriented->end) * .5;
    if (surfaceJet(op, surface, uv.X(), uv.Y()).n.Dot(oriented->reference.jet(op, .5).n) < 0)
      orientation = TopAbs_REVERSED;
  }
  const auto made = surfacing::makeFace(op, surface, tolerance, orientation);
  surfacing::measure(op, made, held, tolerance, normalFlow);
  return stored(k, made);
}
} // namespace BRepExt
