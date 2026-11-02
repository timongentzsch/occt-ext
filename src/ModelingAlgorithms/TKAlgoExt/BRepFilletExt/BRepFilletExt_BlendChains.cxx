// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GC_MakeCircle.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <limits>

namespace BRepExt
{
namespace
{
struct Side
{
  TopoDS_Edge edge;
  int         neighbour;
  int         axis;
  double      at;
  bool        tangent;
  bool        g2;
};

struct Face
{
  int               kind = 0;
  std::vector<Side> sides;
  std::vector<int>  springs;
};

// Circular sections may be delivered as approximated splines. Conic and G2 sections
// retain their quadratic or quintic basis; G2 is not mathematically a conic.
// A quintic may match nonzero support curvature, as on a spherical rim.
bool section(Handle(Geom_Curve) curve, double first, double last, double tolerance, bool g2)
{
  if (!std::isfinite(first) || !std::isfinite(last) || last <= first)
    return false;
  GeomAdaptor_Curve adaptor(curve, first, last);
  const auto        type = adaptor.GetType();
  if (type == GeomAbs_Circle || type == GeomAbs_Ellipse || type == GeomAbs_Parabola
      || type == GeomAbs_Hyperbola)
    return true;
  auto basis = curve;
  while (auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(basis))
    basis = trimmed->BasisCurve();
  int  degree     = 0;
  bool singleSpan = false;
  if (auto spline = Handle(Geom_BSplineCurve)::DownCast(basis))
  {
    degree     = spline->Degree();
    singleSpan = spline->NbKnots() == 2;
  }
  if (auto bezier = Handle(Geom_BezierCurve)::DownCast(basis))
  {
    degree     = bezier->Degree();
    singleSpan = true;
  }
  const auto a = curve->Value(first), b = curve->Value((first + last) / 2), c = curve->Value(last);
  auto       normal = gp_Vec(a, b).Crossed(gp_Vec(a, c));
  if (normal.Magnitude() < tolerance * tolerance)
    return false;
  normal.Normalize();
  const double slack = std::max(tolerance * 10, a.Distance(c) * 1e-5);
  for (int i = 0; i <= 16; ++i)
    if (std::abs(gp_Vec(a, curve->Value(first + (last - first) * i / 16)).Dot(normal)) > slack)
      return false;
  if (degree == 2 && singleSpan)
    return true;
  if (degree == 5 && g2)
    return true;
  GC_MakeCircle circle(a, b, c);
  if (!circle.IsDone())
    return false;
  const auto   round  = circle.Value();
  const double radius = round->Radius();
  for (int i = 0; i <= 16; ++i)
    if (std::abs(curve->Value(first + (last - first) * i / 16).Distance(round->Location()) - radius)
        > slack)
      return false;
  return true;
}

struct Contact
{
  bool tangent;
  bool g2;
};

Contact contact(const TopoDS_Edge& edge,
                const TopoDS_Face& a,
                const TopoDS_Face& b,
                double             tolerance)
{
  BRepAdaptor_Curve curve(edge);
  bool              g2 = true;
  for (int i = 0; i < 9; ++i)
  {
    const double t =
      curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * (i + .5) / 9;
    const auto one   = surfaceJet("blendTopology", edge, a, t),
               other = surfaceJet("blendTopology", edge, b, t);
    if (one.p.Distance(other.p) > std::max(tolerance, BRep_Tool::Tolerance(edge)) * 5
        || one.n.Dot(other.n) < std::cos(M_PI / 360))
      return {false, false};
    gp_Pnt p;
    gp_Vec along;
    curve.D1(t, p, along);
    const auto difference = bendDifference(one, other, along);
    g2 &= difference.first <= 1e-5 || difference.second <= 1e-3;
  }
  return {true, g2};
}
} // namespace

// Rows are face index, kind (strip or corner), count, and cross-edge neighbour indices.
// Indices follow the body's face map, so this read has no dependence on operation history.
std::vector<double> blendTopology(OcctKernel& kernel, uint32_t id, double tolerance)
{
  const auto& shape = get(kernel, id);
  ShapeSet    faces;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  Adjacency adjacent;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacent);
  std::vector<Face> info(faces.Extent());
  for (int index = 1; index <= faces.Extent(); ++index)
  {
    const auto          face  = TopoDS::Face(faces(index));
    auto&               entry = info[index - 1];
    BRepAdaptor_Surface surface(face);
    if (surface.GetType() == GeomAbs_Plane
        || GeomLib_IsPlanarSurface(BRep_Tool::Surface(face), tolerance).IsPlanar())
      continue;
    ShapeSet edges;
    TopExp::MapShapes(face, TopAbs_EDGE, edges);
    for (const auto& e : edges)
    {
      const auto edge = TopoDS::Edge(e);
      if (BRep_Tool::Degenerated(edge))
        continue;
      int neighbour = 0;
      for (const auto& other : adjacent.FindFromKey(edge))
        if (!other.IsSame(face))
        {
          if (neighbour)
          {
            neighbour = 0;
            break;
          }
          neighbour = faces.FindIndex(other);
        }
      if (!neighbour)
        continue;
      double     first, last;
      const auto pc = BRep_Tool::CurveOnSurface(edge, face, first, last);
      if (pc.IsNull())
        continue;
      const auto   a = pc->Value(first), b = pc->Value(last), mid = pc->Value((first + last) / 2);
      const double uv   = std::max(a.Distance(b), 1.0) * 1e-7;
      int          axis = -1;
      if (std::abs(a.X() - b.X()) < uv && std::abs(a.X() - mid.X()) < uv)
        axis = 0;
      else if (std::abs(a.Y() - b.Y()) < uv && std::abs(a.Y() - mid.Y()) < uv)
        axis = 1;
      const auto joint = contact(edge, face, TopoDS::Face(faces(neighbour)), tolerance);
      entry.sides.push_back(
        {edge, neighbour, axis, axis == 0 ? mid.X() : mid.Y(), joint.tangent, joint.g2});
    }
    // Two opposed longitudinal boundaries leave the two supports tangentially.
    // The remaining boundaries are the strip's cross edges, even when sharp.
    if (surface.GetType() == GeomAbs_Sphere)
      continue;
    double shortest = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < entry.sides.size(); ++i)
      for (size_t j = i + 1; j < entry.sides.size(); ++j)
      {
        const auto& a = entry.sides[i];
        const auto& b = entry.sides[j];
        if (!a.tangent || !b.tangent || a.neighbour == b.neighbour || a.axis < 0 || a.axis != b.axis
            || std::abs(a.at - b.at) < 1e-9)
          continue;
        const auto geometry = BRep_Tool::Surface(face);
        bool       valid    = true;
        double     length   = 0;
        for (double fraction : {.2, .5, .8})
        {
          const double along =
            a.axis == 0 ? surface.FirstVParameter()
                            + fraction * (surface.LastVParameter() - surface.FirstVParameter())
                        : surface.FirstUParameter()
                            + fraction * (surface.LastUParameter() - surface.FirstUParameter());
          const auto cross = a.axis == 0 ? geometry->VIso(along) : geometry->UIso(along);
          const auto start = cross->Value(a.at), end = cross->Value(b.at);
          length += start.Distance(end);
          // Opposed supports describe a full wall or a half-round, not a strip.
          const auto jetA = surfaceJet("blendTopology",
                                       geometry,
                                       a.axis == 0 ? a.at : along,
                                       a.axis == 0 ? along : a.at);
          const auto jetB = surfaceJet("blendTopology",
                                       geometry,
                                       a.axis == 0 ? b.at : along,
                                       a.axis == 0 ? along : b.at);
          valid &= jetA.n.Dot(jetB.n) > -1 + 1e-6;
          valid &=
            section(cross, std::min(a.at, b.at), std::max(a.at, b.at), tolerance, a.g2 && b.g2);
        }
        if (valid && length < shortest)
        {
          shortest      = length;
          entry.kind    = 1;
          entry.springs = {int(i), int(j)};
        }
      }
  }
  // A corner patch has three or more strip contacts, rather than two spring edges.
  for (auto& entry : info)
  {
    std::vector<int> neighbours;
    for (const auto& side : entry.sides)
      if (side.tangent && info[side.neighbour - 1].kind)
        neighbours.push_back(side.neighbour);
    std::sort(neighbours.begin(), neighbours.end());
    neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
    if (neighbours.size() >= 3)
    {
      entry.kind = 2;
      entry.springs.clear();
    }
  }
  std::vector<double> out;
  for (size_t i = 0; i < info.size(); ++i)
  {
    const auto& entry = info[i];
    if (!entry.kind)
      continue;
    std::vector<int> cross;
    for (size_t j = 0; j < entry.sides.size(); ++j)
    {
      const auto& side = entry.sides[j];
      if (std::find(entry.springs.begin(), entry.springs.end(), int(j)) == entry.springs.end()
          && info[side.neighbour - 1].kind)
        cross.push_back(side.neighbour - 1);
    }
    std::sort(cross.begin(), cross.end());
    cross.erase(std::unique(cross.begin(), cross.end()), cross.end());
    out.insert(out.end(), {double(i), double(entry.kind), double(cross.size())});
    for (int next : cross)
      out.push_back(next);
  }
  return out;
}
} // namespace BRepExt
