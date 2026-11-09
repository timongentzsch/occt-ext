// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepFilletExt_BlendChecks.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <Geom_BSplineSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Cylinder.hxx>
#include <array>
#include <optional>

namespace BRepExt
{

bool sameBlendSurface(const TopoDS_Face& a, const TopoDS_Face& b)
{
  if (a.IsSame(b))
    return true;
  TopLoc_Location x, y;
  const auto      sa = BRep_Tool::Surface(a, x), sb = BRep_Tool::Surface(b, y);
  if (sa == sb && x.IsEqual(y))
    return true;
  const BRepAdaptor_Surface aa(a, false), bb(b, false);
  if (aa.GetType() != bb.GetType())
    return false;
  if (aa.GetType() == GeomAbs_Plane)
    return aa.Plane().Axis().Direction().IsParallel(bb.Plane().Axis().Direction(), tolerance::angle)
           && aa.Plane().Distance(bb.Plane().Location()) < tolerance::length;
  if (aa.GetType() == GeomAbs_Cylinder)
  {
    const auto p = aa.Cylinder(), q = bb.Cylinder();
    return std::abs(p.Radius() - q.Radius()) < tolerance::length
           && p.Axis().Direction().IsParallel(q.Axis().Direction(), tolerance::angle)
           && gp_Lin(p.Axis()).Distance(q.Location()) < tolerance::length;
  }
  return false;
}

bool tangentBlendEdge(const TopoDS_Shape& shape, const TopoDS_Edge& edge)
{
  if (BRep_Tool::Degenerated(edge))
    return true;
  Adjacency adjacent;
  TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacent);
  if (adjacent.Contains(edge) && adjacent.FindFromKey(edge).Extent() == 1
      && BRep_Tool::IsClosed(edge, TopoDS::Face(adjacent.FindFromKey(edge).First())))
    return true;
  if (!adjacent.Contains(edge) || adjacent.FindFromKey(edge).Extent() != 2)
    fail("fillet",
         unsupported_geometry,
         "an edge to blend must lie between two distinct faces",
         {edge},
         FailureKind::InvalidTopology);
  const auto&       pair = adjacent.FindFromKey(edge);
  const BoundaryJet a("fillet", edge, TopoDS::Face(pair.First())),
    b("fillet", edge, TopoDS::Face(pair.Last()));
  BRepAdaptor_Curve curve(edge);
  for (int i = 0; i < 9; ++i)
  {
    const double t =
      curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * (i + .5) / 9;
    const auto x = a.at(t), y = b.at(t);
    if (x.n.Dot(y.n) <= 0 || x.n.Crossed(y.n).Magnitude() > tolerance::coarseAngle)
      return false;
  }
  return true;
}

BlendBoundaryRoles blendBoundaryRoles(const TopoDS_Shape&              before,
                                      const std::vector<TopoDS_Shape>& selected,
                                      const TopoDS_Shape&              result,
                                      const BooleanResult&             history,
                                      const std::vector<TopoDS_Shape>& generated)
{
  ShapeSet supports;
  for (TopExp_Explorer f(before, TopAbs_FACE); f.More(); f.Next())
    for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
      for (const auto& pick : selected)
        if (e.Current().IsSame(pick))
          supports.Add(f.Current());
  BlendBoundaryRoles roles;
  for (TopExp_Explorer f(result, TopAbs_FACE); f.More(); f.Next())
  {
    bool fresh = false;
    for (const auto& made : generated)
      fresh |= f.Current().IsSame(made);
    if (fresh)
    {
      roles.blends.Add(f.Current());
      continue;
    }
    bool inherited = false, support = false;
    for (TopExp_Explorer old(before, TopAbs_FACE); old.More(); old.Next())
    {
      bool match = sameBlendSurface(TopoDS::Face(f.Current()), TopoDS::Face(old.Current()));
      if (const auto* descendants = history.faces.Seek(old.Current()))
        for (const auto& descendant : *descendants)
          match |= f.Current().IsSame(descendant);
      if (match)
      {
        inherited = true;
        support |= supports.Contains(old.Current());
      }
    }
    if (!inherited)
      roles.blends.Add(f.Current());
    else
      (support ? roles.supports : roles.stops).Add(f.Current());
  }
  return roles;
}

void certifyBlendBoundaries(const char*               op,
                            const TopoDS_Shape&       result,
                            const BlendBoundaryRoles& roles,
                            double                    scale,
                            double                    angleDegrees,
                            int                       samples)
{
  Adjacency adjacency;
  TopExp::MapShapesAndUniqueAncestors(result, TopAbs_EDGE, TopAbs_FACE, adjacency);
  const double cosine = std::cos(angleDegrees * M_PI / 180);
  for (int i = 1; i <= adjacency.Extent(); ++i)
  {
    const auto  edge  = TopoDS::Edge(adjacency.FindKey(i));
    const auto& pair  = adjacency.FindFromIndex(i);
    bool        blend = false;
    for (const auto& face : pair)
      blend |= roles.blends.Contains(face);
    if (!blend || BRep_Tool::Degenerated(edge))
      continue;
    if (pair.Extent() == 1 && BRep_Tool::IsClosed(edge, TopoDS::Face(pair.First())))
      continue;
    if (pair.Extent() < 2)
      fail(op,
           no_solution,
           "blend boundary has no adjoining face",
           {edge},
           FailureKind::InvalidTopology);
    // OCCT can retain a tangent three-face contact on a valid closed solid.
    // Certify every face pair there, instead of rejecting it by arity or
    // skipping the junction as the older two-face-only checks did.
    std::vector<TopoDS_Face> faces;
    for (const auto& face : pair)
      faces.push_back(TopoDS::Face(face));
    BRepAdaptor_Curve curve(edge);
    for (size_t first = 0; first < faces.size(); ++first)
      for (size_t second = first + 1; second < faces.size(); ++second)
      {
        const auto a = faces[first], b = faces[second];
        bool       stop =
          roles.sharp.Contains(edge) || roles.stops.Contains(a) || roles.stops.Contains(b);
        for (const auto& [stripe, contacts] : roles.contacts)
        {
          if (stripe.IsSame(a) && !roles.blends.Contains(b) && !contacts.Contains(b))
            stop = true;
          if (stripe.IsSame(b) && !roles.blends.Contains(a) && !contacts.Contains(a))
            stop = true;
        }
        // Honor both the builder's approximation contract and declared shape
        // tolerances when comparing jets with the 3D edge;
        // the G1 angle bound is independent and never widened by this tolerance.
        const double      contactTolerance = std::max({roles.constructionTolerance,
                                                       tolerance::sewing,
                                                       1e-5 * scale,
                                                       BRep_Tool::Tolerance(edge),
                                                       BRep_Tool::Tolerance(a),
                                                       BRep_Tool::Tolerance(b)});
        const char*       role             = stop ? "sharp stop"
                                             : roles.blends.Contains(a) && roles.blends.Contains(b)
                                               ? "blend seam"
                                               : "support contact";
        const BoundaryJet ours(op, edge, a), theirs(op, edge, b);
        std::optional<std::array<GeomAPI_ProjectPointOnSurf, 2>> stopProjections;
        for (int j = 0; j < samples; ++j)
        {
          const double t = curve.FirstParameter()
                           + (curve.LastParameter() - curve.FirstParameter()) * (j + .5) / samples;
          const auto x = ours.at(t), y = theirs.at(t);
          const auto p = curve.Value(t);
          if (std::max({x.p.Distance(y.p), x.p.Distance(p), y.p.Distance(p)}) > contactTolerance)
          {
            // Sharp miters need not share the 3D curve's parameterization.
            // Certify the actual curve on both surfaces in that case, at
            // the same position bound. Smooth contacts still require
            // matching jets for their G1 reading.
            bool onStop = stop;
            if (stop)
            {
              if (!stopProjections)
              {
                stopProjections.emplace();
                const TopoDS_Face sides[2] = {a, b};
                for (int side = 0; side < 2; ++side)
                {
                  double u0, u1, v0, v1;
                  BRepTools::UVBounds(sides[side], u0, u1, v0, v1);
                  (*stopProjections)[side]
                    .Init(BRep_Tool::Surface(sides[side]), u0, u1, v0, v1, tolerance::uv);
                }
              }
              for (auto& projection : *stopProjections)
              {
                projection.Perform(p);
                onStop &=
                  projection.NbPoints() > 0 && projection.LowerDistance() <= contactTolerance;
              }
            }
            if (!onStop)
              fail(op,
                   no_solution,
                   std::string(role) + " leaves its shared boundary curve",
                   {edge},
                   FailureKind::SupportContact);
          }
          if (!stop && std::abs(x.n.Dot(y.n)) < cosine)
            fail(op,
                 no_solution,
                 std::string(role) + " is not G1: normal mismatch "
                   + std::to_string(std::acos(std::clamp(std::abs(x.n.Dot(y.n)), 0.0, 1.0)) * 180
                                    / M_PI)
                   + " degrees",
                 {edge},
                 FailureKind::Tangency);
        }
      }
  }
}
} // namespace BRepExt
