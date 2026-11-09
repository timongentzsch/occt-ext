// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Give one face of a body a new support and let its neighbours meet it, as when the control
// points of a face are moved:
//
//   - each corner of the face slides along the edges that leave the face there, to where they
//     cross the new support;
//   - each edge of the face becomes the curve where the new support meets the neighbour across
//     it, between the new corners;
//   - the neighbours keep their supports and are grown or cut to the new edges.
//
// A corner or an edge that still lies on the new support stays as it is, so moving inner
// control points leaves the neighbours alone. Where the new support does not reach the face's
// new boundary it is extended along its own boundary rows (GeomLib's extension), and the face is
// given that much of the extension: the support then has spans the net that was asked for does
// not, which is why a nominal surface is kept.
//
// The topology does not change. An edit that would need it to (a corner that leaves its edge,
// an edge that the neighbour no longer meets, a free or seam edge, a neighbour whose support
// ends short of the new edge) is refused, naming the face or edge.

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepSurfaceExt_Tools.hxx>

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_HArray1.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <cmath>
#include <limits>

namespace BRepExt
{
namespace
{

using ShapeMap = NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>;

occ::handle<Geom_Surface> worldSurface(const TopoDS_Face& face)
{
  TopLoc_Location location;
  auto            surface = BRep_Tool::Surface(face, location);
  if (location.IsIdentity())
    return surface;
  return occ::down_cast<Geom_Surface>(surface->Transformed(location.Transformation()));
}

//! The curve an edge lies on, in world coordinates, and beyond its ends where its kind goes on:
//! a line or a circle whole, a spline as far as it is defined.
occ::handle<Geom_Curve> carrier(const TopoDS_Edge& edge, double& first, double& last)
{
  TopLoc_Location location;
  auto            curve = BRep_Tool::Curve(edge, location, first, last);
  if (curve.IsNull())
    return curve;
  while (auto trimmed = occ::down_cast<Geom_TrimmedCurve>(curve))
    curve = trimmed->BasisCurve();
  if (!location.IsIdentity())
    curve = occ::down_cast<Geom_Curve>(curve->Transformed(location.Transformation()));
  return curve;
}

double distanceTo(const gp_Pnt& point, const occ::handle<Geom_Surface>& surface)
{
  GeomAPI_ProjectPointOnSurf projection(point, surface);
  return projection.NbPoints() ? projection.LowerDistance()
                               : std::numeric_limits<double>::infinity();
}

//! The largest distance from `count` points along a curve's range to `surface`.
double offSurface(const occ::handle<Geom_Curve>&   curve,
                  double                           first,
                  double                           last,
                  const occ::handle<Geom_Surface>& surface,
                  int                              count = 16)
{
  double worst = 0;
  for (int i = 0; i <= count; ++i)
    worst = std::max(worst, distanceTo(curve->Value(first + (last - first) * i / count), surface));
  return worst;
}

//! Where two surfaces meet near `near`, in the plane through it across `along` (a unit
//! vector): Newton's method on both surfaces' parameters, from the nearest points of each.
bool meetAt(const occ::handle<Geom_Surface>& a,
            const occ::handle<Geom_Surface>& b,
            const gp_Pnt&                    near,
            const gp_Vec&                    along,
            double                           tolerance,
            gp_Pnt&                          out)
{
  GeomAPI_ProjectPointOnSurf onA(near, a), onB(near, b);
  if (!onA.NbPoints() || !onB.NbPoints())
    return false;
  double u, v, s, t;
  onA.LowerDistanceParameters(u, v);
  onB.LowerDistanceParameters(s, t);
  for (int step = 0; step < 50; ++step)
  {
    gp_Pnt pa, pb;
    gp_Vec au, av, bs, bt;
    a->D1(u, v, pa, au, av);
    b->D1(s, t, pb, bs, bt);
    const gp_Vec gap(pb, pa);
    const double across = gp_Vec(near, pa).Dot(along);
    if (gap.Magnitude() < 1e-3 * tolerance && std::abs(across) < 1e-3 * tolerance)
    {
      out = pa;
      return true;
    }
    math_Matrix jacobian(1, 4, 1, 4, 0.0);
    math_Vector rhs(1, 4), delta(1, 4);
    for (int r = 1; r <= 3; ++r)
    {
      jacobian(r, 1) = au.Coord(r);
      jacobian(r, 2) = av.Coord(r);
      jacobian(r, 3) = -bs.Coord(r);
      jacobian(r, 4) = -bt.Coord(r);
      rhs(r)         = -gap.Coord(r);
    }
    jacobian(4, 1) = au.Dot(along);
    jacobian(4, 2) = av.Dot(along);
    rhs(4)         = -across;
    math_Gauss solve(jacobian);
    if (!solve.IsDone())
      return false;
    solve.Solve(rhs, delta);
    u += delta(1), v += delta(2), s += delta(3), t += delta(4);
  }
  return false;
}

//! The curve where `a` and `b` meet, from `start` to `end` (once round, when `closed`), found
//! across the old edge's curve point by point and interpolated. Points are added between two
//! where the curve strays from either surface by more than half the tolerance.
occ::handle<Geom_Curve> meetingCurve(const char*                      op,
                                     const occ::handle<Geom_Curve>&   old,
                                     double                           first,
                                     double                           last,
                                     const occ::handle<Geom_Surface>& a,
                                     const occ::handle<Geom_Surface>& b,
                                     const gp_Pnt&                    start,
                                     const gp_Pnt&                    end,
                                     bool                             closed,
                                     double                           tolerance,
                                     const std::vector<TopoDS_Shape>& at)
{
  struct Site
  {
    double t;
    gp_Pnt point;
  };

  auto solve = [&](double t) {
    gp_Pnt near;
    gp_Vec along;
    old->D1(t, near, along);
    gp_Pnt out;
    if (along.Magnitude() <= gp::Resolution()
        || !meetAt(a, b, near, along.Normalized(), tolerance, out))
      fail(op,
           no_solution,
           "the new surface does not meet the face beside this edge, or meets it tangentially",
           at);
    return Site{t, out};
  };
  std::vector<Site> sites{{first, start}};
  const int         count = closed ? 32 : 16;
  for (int i = 1; i < count; ++i)
    sites.push_back(solve(first + (last - first) * i / count));
  if (!closed)
    sites.push_back({last, end});
  for (;;)
  {
    occ::handle<NCollection_HArray1<gp_Pnt>> points =
      new NCollection_HArray1<gp_Pnt>(1, int(sites.size()));
    for (size_t i = 0; i < sites.size(); ++i)
      points->SetValue(int(i) + 1, sites[i].point);
    GeomAPI_Interpolate interpolate(points, closed, Precision::Confusion());
    interpolate.Perform();
    if (!interpolate.IsDone())
      fail(op, no_solution, "could not run a curve through where the surfaces meet", at);
    const occ::handle<Geom_BSplineCurve> curve = interpolate.Curve();
    // Between each two sites: the curve's distance from both surfaces.
    std::vector<Site> refined;
    bool              good  = true;
    const int         spans = closed ? int(sites.size()) : int(sites.size()) - 1;
    for (int i = 0; i < spans; ++i)
    {
      refined.push_back(sites[i]);
      const double low = curve->Knot(i + 1), high = curve->Knot(i + 2);
      const gp_Pnt middle = curve->Value((low + high) / 2);
      if (std::max(distanceTo(middle, a), distanceTo(middle, b)) > tolerance / 2)
      {
        const double next = i + 1 < int(sites.size()) ? sites[i + 1].t : last;
        refined.push_back(solve((sites[i].t + next) / 2));
        good = false;
      }
    }
    if (!closed)
      refined.push_back(sites.back());
    if (good)
      return Handle(Geom_Curve)(curve);
    if (refined.size() > 4096)
      fail(op, no_solution, "the curve where the surfaces meet is too fine to follow", at);
    sites = std::move(refined);
  }
}

struct Edge
{
  TopoDS_Edge             original; // in the copy, forward
  TopoDS_Edge             made;     // its replacement, forward
  occ::handle<Geom_Curve> curve;    // world, of `made`
  double                  first = 0, last = 0;
};

} // namespace

EvolutionData replaceFaceSupport(OcctKernel&                k,
                                 uint32_t                   bodyId,
                                 uint32_t                   faceId,
                                 const std::vector<double>& spline,
                                 double                     tolerance,
                                 const std::vector<int>&    hashes,
                                 int                        bound)
{
  const char* op = "replaceFaceSupport";
  if (!std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "positive tolerance required");
  const auto& input = get(k, bodyId);
  ShapeSet    members;
  TopExp::MapShapes(input, TopAbs_FACE, members);
  const auto& picked = get(k, faceId);
  if (picked.ShapeType() != TopAbs_FACE || !members.Contains(picked))
    fail(op, invalid_input, "select a face of the body");
  const auto wanted = splineFromData(op, spline);
  if (wanted->IsUPeriodic() || wanted->IsVPeriodic())
    fail(op, unsupported_geometry, "periodic nets are unsupported");
  FaceHistory history(op, hashes, bound);

  // Everything below works on a copy: tolerances are raised in place.
  BRepBuilderAPI_Copy copy(input, true, false);
  const auto          body = copy.Shape();
  const auto face = TopoDS::Face(copy.ModifiedShape(members.FindKey(members.FindIndex(picked))));
  // Failures name the input's own sub-shapes, which core knows.
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> origin;
  for (const auto type : {TopAbs_EDGE, TopAbs_VERTEX, TopAbs_FACE})
    for (TopExp_Explorer it(input, type); it.More(); it.Next())
      if (!origin.IsBound(copy.ModifiedShape(it.Current())))
        origin.Bind(copy.ModifiedShape(it.Current()), it.Current());
  auto at = [&](const TopoDS_Shape& shape) {
    return std::vector<TopoDS_Shape>{origin.IsBound(shape) ? origin.Find(shape) : shape};
  };

  Adjacency edgeFaces, vertexEdges;
  TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, edgeFaces);
  TopExp::MapShapesAndAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, vertexEdges);
  ShapeSet faceEdges, faceVertices;
  TopExp::MapShapes(face, TopAbs_EDGE, faceEdges);
  TopExp::MapShapes(face, TopAbs_VERTEX, faceVertices);

  // The one neighbour across each edge.
  NCollection_DataMap<TopoDS_Shape, TopoDS_Face, TopTools_ShapeMapHasher> across;
  for (const auto& shape : faceEdges)
  {
    const auto edge = TopoDS::Edge(shape);
    if (BRep_Tool::Degenerated(edge) || BRep_Tool::IsClosed(edge, face))
      fail(op,
           unsupported_geometry,
           "faces with seams or degenerate edges are unsupported",
           at(edge));
    std::vector<TopoDS_Face> others;
    for (const auto& owner : edgeFaces.FindFromKey(edge))
      if (!owner.IsSame(face))
        others.push_back(TopoDS::Face(owner));
    if (others.size() != 1)
      fail(op,
           unsupported_geometry,
           others.empty() ? "an edge with no neighbour is unsupported"
                          : "an edge with several neighbours is unsupported",
           at(edge));
    across.Bind(edge, others.front());
  }

  // As far as the new boundary can be from the old one: the face's own size.
  Bnd_Box box;
  BRepBndLib::Add(face, box);
  const double                     reach = std::sqrt(box.SquareExtent());
  occ::handle<Geom_BoundedSurface> grown = occ::down_cast<Geom_BSplineSurface>(wanted->Copy());
  for (const bool alongU : {true, false})
    for (const bool after : {true, false})
      GeomLib::ExtendSurfByLength(grown,
                                  reach,
                                  std::min(wanted->UDegree(), wanted->VDegree()) >= 2 ? 2 : 1,
                                  alongU,
                                  after);
  const auto extended = occ::down_cast<Geom_BSplineSurface>(grown);
  if (extended.IsNull())
    fail(op, no_solution, "could not extend the new support");

  // Corners: kept where they are on the new support, else slid along the edges that leave
  // the face there.
  ShapeMap     corners;
  BRep_Builder builder;
  for (const auto& shape : faceVertices)
  {
    const auto vertex = TopoDS::Vertex(shape);
    const auto point  = BRep_Tool::Pnt(vertex);
    if (distanceTo(point, extended) <= tolerance)
    {
      corners.Bind(vertex, vertex);
      continue;
    }
    std::vector<gp_Pnt> found;
    for (const auto& incident : vertexEdges.FindFromKey(vertex))
    {
      if (faceEdges.Contains(incident))
        continue;
      const auto side = TopoDS::Edge(incident);
      if (BRep_Tool::Degenerated(side))
        continue;
      double     first, last;
      const auto curve = carrier(side, first, last);
      if (curve.IsNull())
        fail(op, no_solution, "an edge has no curve", at(side));
      GeomAPI_IntCS crossing(curve, extended);
      if (!crossing.IsDone() || crossing.NbPoints() == 0)
        fail(op, no_solution, "an edge leaving the face does not reach the new surface", at(side));
      int best = 1;
      for (int i = 2; i <= crossing.NbPoints(); ++i)
        if (crossing.Point(i).Distance(point) < crossing.Point(best).Distance(point))
          best = i;
      found.push_back(crossing.Point(best));
    }
    if (found.empty())
      fail(op,
           unsupported_geometry,
           "a corner with no edge to slide along is unsupported",
           at(vertex));
    for (const auto& other : found)
      if (other.Distance(found.front()) > tolerance)
        fail(op,
             unsupported_geometry,
             "the edges at a corner no longer meet on the new surface",
             at(vertex));
    TopoDS_Vertex moved;
    builder.MakeVertex(moved, found.front(), BRep_Tool::Tolerance(vertex));
    corners.Bind(vertex, moved);
  }
  auto cornerOf = [&](const TopoDS_Vertex& vertex) { return TopoDS::Vertex(corners.Find(vertex)); };

  // Where a curve passes `point`: its parameter.
  auto parameterOn =
    [&](const occ::handle<Geom_Curve>& curve, const gp_Pnt& point, double& parameter) {
      GeomAPI_ProjectPointOnCurve projection(point, curve);
      if (projection.NbPoints() == 0)
        return std::numeric_limits<double>::infinity();
      parameter = projection.LowerDistanceParameter();
      return projection.LowerDistance();
    };
  // A new edge on `curve` between two corners, as `original` ran. An edge that leaves the
  // face (`turnable` false) keeps its curve's direction: a corner that slid past its far end
  // would turn it round.
  auto edgeOn = [&](const TopoDS_Edge&      original,
                    occ::handle<Geom_Curve> curve,
                    const TopoDS_Vertex&    start,
                    const TopoDS_Vertex&    end,
                    bool                    turnable) {
    double a, b;
    if (parameterOn(curve, BRep_Tool::Pnt(start), a) > 10 * tolerance
        || parameterOn(curve, BRep_Tool::Pnt(end), b) > 10 * tolerance)
      fail(op, no_solution, "the new edge does not pass through its corners", at(original));
    const bool closed = start.IsSame(end);
    if (closed)
    {
      // A closed intersection comes back as a spline that ends where it starts; as a
      // periodic one its seam can be where the corner is.
      if (auto spline = occ::down_cast<Geom_BSplineCurve>(curve);
          !curve->IsPeriodic() && !spline.IsNull() && spline->IsClosed())
      {
        spline = occ::down_cast<Geom_BSplineCurve>(spline->Copy());
        try
        {
          spline->SetPeriodic();
        }
        catch (const Standard_Failure&)
        {
          fail(op,
               unsupported_geometry,
               "a closed edge on an open curve is unsupported",
               at(original));
        }
        curve = spline;
        parameterOn(curve, BRep_Tool::Pnt(start), a);
      }
      if (!curve->IsPeriodic())
        fail(op,
             unsupported_geometry,
             "a closed edge on an open curve is unsupported",
             at(original));
      b = a + curve->Period();
    }
    // The new edge runs the way the old one did, which the old midpoint tangent says.
    double     first, last;
    const auto old = carrier(original, first, last);
    gp_Pnt     middle;
    gp_Vec     was, now;
    old->D1((first + last) / 2, middle, was);
    double near;
    parameterOn(curve, middle, near);
    gp_Pnt on;
    curve->D1(near, on, now);
    if (!closed && b < a)
    {
      if (curve->IsPeriodic())
        b += curve->Period() * std::ceil((a - b) / curve->Period());
      else if (!turnable)
        fail(op,
             unsupported_geometry,
             "a corner would slide past the end of its edge",
             at(original));
      else
      {
        curve = curve->Reversed();
        a     = curve->ReversedParameter(a);
        b     = curve->ReversedParameter(b);
      }
    }
    if (closed && was.Dot(now) < 0)
    {
      curve = curve->Reversed();
      a     = curve->ReversedParameter(a);
      b     = a + curve->Period();
    }
    if (b - a <= Precision::PConfusion())
      fail(op, unsupported_geometry, "an edge would vanish", at(original));
    TopoDS_Edge made;
    builder.MakeEdge(made, curve, BRep_Tool::Tolerance(original));
    builder.Add(made, start.Oriented(TopAbs_FORWARD));
    builder.Add(made, end.Oriented(TopAbs_REVERSED));
    builder.Range(made, a, b);
    return Edge{original, made, curve, a, b};
  };

  // The face's edges.
  std::vector<Edge> edges;
  for (const auto& shape : faceEdges)
  {
    const auto    edge = TopoDS::Edge(shape.Oriented(TopAbs_FORWARD));
    TopoDS_Vertex v0, v1;
    TopExp::Vertices(edge, v0, v1);
    double     first, last;
    const auto curve = carrier(edge, first, last);
    const auto start = cornerOf(v0), end = cornerOf(v1);
    if (start.IsSame(v0) && end.IsSame(v1) && offSurface(curve, first, last, extended) <= tolerance)
    {
      edges.push_back(edgeOn(edge, curve, start, end, false));
      continue;
    }
    const auto curve2 = meetingCurve(op,
                                     curve,
                                     first,
                                     last,
                                     extended,
                                     worldSurface(across.Find(edge)),
                                     BRep_Tool::Pnt(start),
                                     BRep_Tool::Pnt(end),
                                     start.IsSame(end),
                                     tolerance,
                                     at(edge));
    edges.push_back(edgeOn(edge, curve2, start, end, true));
  }

  // The edges that leave the face at a moved corner run to it now.
  std::vector<Edge> sides;
  ShapeSet          touched;
  for (const auto& shape : faceVertices)
  {
    const auto vertex = TopoDS::Vertex(shape);
    if (cornerOf(vertex).IsSame(vertex))
      continue;
    for (const auto& incident : vertexEdges.FindFromKey(vertex))
      if (!faceEdges.Contains(incident) && !BRep_Tool::Degenerated(TopoDS::Edge(incident)))
        touched.Add(incident.Oriented(TopAbs_FORWARD));
  }
  for (const auto& shape : touched)
  {
    const auto    side = TopoDS::Edge(shape);
    TopoDS_Vertex v0, v1;
    TopExp::Vertices(side, v0, v1);
    const auto start = corners.IsBound(v0) ? cornerOf(v0) : v0,
               end   = corners.IsBound(v1) ? cornerOf(v1) : v1;
    double first, last;
    sides.push_back(edgeOn(side, carrier(side, first, last), start, end, false));
  }

  // The support the face gets: the new net where it reaches the new boundary, else the net and
  // as much of its extension as the boundary needs, up to the boundary, so that a face the net
  // has shrunk inside is again a whole face of its support.
  double u0, u1, v0, v1;
  wanted->Bounds(u0, u1, v0, v1);
  double lowU = u0, highU = u1, lowV = v0, highV = v1;
  bool   beyond = false;
  for (const auto& edge : edges)
    for (int i = 0; i <= 32; ++i)
    {
      GeomAPI_ProjectPointOnSurf projection(
        edge.curve->Value(edge.first + (edge.last - edge.first) * i / 32),
        extended);
      if (!projection.NbPoints())
        fail(op, no_solution, "an edge left the extended surface", at(edge.original));
      double u, v;
      projection.LowerDistanceParameters(u, v);
      const double marginU = 1e-9 * (u1 - u0), marginV = 1e-9 * (v1 - v0);
      if (u < u0 - marginU || u > u1 + marginU || v < v0 - marginV || v > v1 + marginV)
        beyond = true;
      lowU = std::min(lowU, u), highU = std::max(highU, u), lowV = std::min(lowV, v),
      highV = std::max(highV, v);
    }
  occ::handle<Geom_BSplineSurface> support = wanted;
  if (beyond)
  {
    support = occ::down_cast<Geom_BSplineSurface>(extended->Copy());
    double e0, e1, f0, f1;
    support->Bounds(e0, e1, f0, f1);
    support->Segment(std::max(e0, lowU),
                     std::min(e1, highU),
                     std::max(f0, lowV),
                     std::min(f1, highV));
  }

  // The new face, its edges' p-curves, and the neighbours' p-curves of every new edge.
  TopoDS_Face made;
  builder.MakeFace(made, support, TopLoc_Location(), BRep_Tool::Tolerance(face));
  // An edge's p-curve on `surface`: its points projected, carried across the surface's period
  // so they run on without a jump, and interpolated at the edge's own parameters, so curve and
  // p-curve agree point for point. On a face that the edge's original lay on, it starts in
  // the same period as the original's p-curve did.
  auto laid = [&](const Edge& edge, const occ::handle<Geom_Surface>& surface, const gp_Pnt2d* was) {
    const double uPeriod = surface->IsUPeriodic() ? surface->UPeriod() : 0,
                 vPeriod = surface->IsVPeriodic() ? surface->VPeriod() : 0;
    auto wrap            = [](double value, double toward, double period) {
      return period > 0 ? nearPeriod(value, toward, period) : value;
    };
    for (int count = 32;; count *= 2)
    {
      occ::handle<NCollection_HArray1<gp_Pnt2d>> points =
        new NCollection_HArray1<gp_Pnt2d>(1, count + 1);
      occ::handle<NCollection_HArray1<double>> parameters =
        new NCollection_HArray1<double>(1, count + 1);
      for (int i = 0; i <= count; ++i)
      {
        const double               t = edge.first + (edge.last - edge.first) * i / count;
        GeomAPI_ProjectPointOnSurf projection(edge.curve->Value(t), surface);
        if (!projection.NbPoints())
          fail(op, no_solution, "could not lay an edge on its face", at(edge.original));
        double u, v;
        projection.LowerDistanceParameters(u, v);
        const gp_Pnt2d toward = i > 0 ? points->Value(i) : was ? *was : gp_Pnt2d(u, v);
        points->SetValue(i + 1,
                         gp_Pnt2d(wrap(u, toward.X(), uPeriod), wrap(v, toward.Y(), vPeriod)));
        parameters->SetValue(i + 1, t);
      }
      Geom2dAPI_Interpolate interpolate(points, parameters, false, Precision::PConfusion());
      interpolate.Perform();
      if (!interpolate.IsDone())
        fail(op, no_solution, "could not lay an edge on its face", at(edge.original));
      const occ::handle<Geom2d_Curve> curve2d(interpolate.Curve());
      double                          worst = 0;
      for (int i = 0; i < count; ++i)
      {
        const double t  = edge.first + (edge.last - edge.first) * (i + 0.5) / count;
        const auto   uv = curve2d->Value(t);
        worst = std::max(worst, surface->Value(uv.X(), uv.Y()).Distance(edge.curve->Value(t)));
      }
      if (worst <= tolerance / 2)
        return curve2d;
      if (count >= 4096)
        fail(op,
             no_solution,
             "could not lay an edge on its face closely enough",
             at(edge.original));
    }
  };
  auto pcurve = [&](Edge& edge, const TopoDS_Face& on, const occ::handle<Geom_Surface>& surface) {
    const double within = std::max(BRep_Tool::Tolerance(edge.made), tolerance);
    // Read as the face is stored: a reversed face would swap a seam's two p-curves.
    const auto     stored = TopoDS::Face(on.Oriented(TopAbs_FORWARD));
    double         f0, l0;
    const auto     forward = on.IsSame(made) ? occ::handle<Geom2d_Curve>()
                                             : BRep_Tool::CurveOnSurface(edge.original, stored, f0, l0);
    const gp_Pnt2d start   = forward.IsNull() ? gp_Pnt2d() : forward->Value(f0);
    const auto     curve2d = laid(edge, surface, forward.IsNull() ? nullptr : &start);
    if (forward.IsNull() || !BRep_Tool::IsClosed(edge.original, on))
    {
      builder.UpdateEdge(edge.made, curve2d, on, within);
      return;
    }
    // A seam has a p-curve on each side of the face's period, as the old one had.
    double     f1, l1;
    const auto reversed =
      BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.original.Reversed()), stored, f1, l1);
    const auto across =
      occ::down_cast<Geom2d_Curve>(curve2d->Translated(gp_Vec2d(start, reversed->Value(f1))));
    builder.UpdateEdge(edge.made, curve2d, across, on, within);
  };
  ShapeMap replaced;
  for (auto& edge : edges)
  {
    pcurve(edge, made, support);
    pcurve(edge, across.Find(edge.original), worldSurface(across.Find(edge.original)));
    replaced.Bind(edge.original, edge.made);
  }
  for (auto& side : sides)
  {
    for (const auto& owner : edgeFaces.FindFromKey(side.original))
    {
      pcurve(side, TopoDS::Face(owner), worldSurface(TopoDS::Face(owner)));
    }
    replaced.Bind(side.original, side.made);
  }
  for (TopExp_Explorer wires(face.Oriented(TopAbs_FORWARD), TopAbs_WIRE); wires.More();
       wires.Next())
  {
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for (TopoDS_Iterator it(wires.Current()); it.More(); it.Next())
      builder.Add(
        wire,
        replaced.Find(it.Value().Oriented(TopAbs_FORWARD)).Oriented(it.Value().Orientation()));
    builder.Add(made, wire.Oriented(wires.Current().Orientation()));
  }
  made.Orientation(face.Orientation());
  for (auto& edge : edges)
    BRepLib::SameParameter(edge.made, tolerance);
  for (auto& side : sides)
    BRepLib::SameParameter(side.made, tolerance);

  BRepTools_ReShape reshape;
  reshape.Replace(face, made);
  for (const auto& edge : edges)
    reshape.Replace(edge.original, edge.made);
  for (const auto& side : sides)
    reshape.Replace(side.original, side.made);
  const auto result = reshape.Apply(body);
  BRepLib::UpdateTolerances(result, false);
  if (!BRepCheck_Analyzer(result).IsValid())
    fail(op,
         no_solution,
         "the faces beside it could not be grown or cut to the new surface",
         at(face));
  for (const auto& each : members)
    if (history.tracks(each))
    {
      NCollection_List<TopoDS_Shape> successor;
      successor.Append(reshape.Value(copy.ModifiedShape(each)));
      history.modified(each, successor);
    }
  return history.finish(store(k, result));
}

} // namespace BRepExt
