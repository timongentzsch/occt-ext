// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <BRepFilletExt_BlendChecks.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <GeomAlgoExt_SurfaceBounds.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BndLib_AddSurface.hxx>
#include <Bnd_Box.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
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

// Two incident stripes terminate in a construction-defined bisector miter.
// The remaining source edge carries the sharp corner through the joint.
void markMiterBoundaries(const TopoDS_Shape&              before,
                         const std::vector<TopoDS_Shape>& selected,
                         const TopoDS_Shape&              result,
                         BlendBoundaryRoles&              roles)
{
  Adjacency vertices, seams;
  TopExp::MapShapesAndUniqueAncestors(before, TopAbs_VERTEX, TopAbs_EDGE, vertices);
  TopExp::MapShapesAndUniqueAncestors(result, TopAbs_EDGE, TopAbs_FACE, seams);
  for (int i = 1; i <= vertices.Extent(); ++i)
  {
    const auto               vertex = TopoDS::Vertex(vertices.FindKey(i));
    const auto               point  = BRep_Tool::Pnt(vertex);
    std::vector<TopoDS_Edge> incident;
    const auto&              edges = vertices.FindFromIndex(i);
    for (const auto& edge : edges)
      for (const auto& pick : selected)
        if (edge.IsSame(pick))
          incident.push_back(TopoDS::Edge(edge));
    if (incident.size() != 2 || edges.Extent() != 3)
      continue;
    gp_Vec direction[2];
    bool   straight = true;
    for (int side = 0; side < 2; ++side)
    {
      BRepAdaptor_Curve curve(incident[side]);
      straight &= curve.GetType() == GeomAbs_Line;
      const auto a = curve.Value(curve.FirstParameter()), b = curve.Value(curve.LastParameter());
      direction[side] = point.Distance(a) < point.Distance(b) ? gp_Vec(a, b) : gp_Vec(b, a);
      direction[side].Normalize();
    }
    if (!straight || direction[0].Crossed(direction[1]).Magnitude() < tolerance::coarseAngle)
      continue;
    const gp_Pln bisector(point, gp_Dir(direction[0] - direction[1]));
    for (int j = 1; j <= seams.Extent(); ++j)
    {
      const auto& pair = seams.FindFromIndex(j);
      if (pair.Extent() != 2 || !roles.blends.Contains(pair.First())
          || !roles.blends.Contains(pair.Last()))
        continue;
      const auto edge = TopoDS::Edge(seams.FindKey(j));
      if (BRep_Tool::Degenerated(edge))
        continue;
      BRepAdaptor_Curve curve(edge);
      bool              on = true;
      for (int k = 0; k <= 8; ++k)
        on &= bisector.Distance(
                curve.Value(curve.FirstParameter()
                            + (curve.LastParameter() - curve.FirstParameter()) * (k + .5) / 9))
              < tolerance::sewing;
      if (on)
        roles.sharp.Add(edge);
    }
  }
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

//! A face the blend made that runs into a face behind its supports is no blend.
void checkBlendClear(const TopoDS_Shape&              before,
                     const TopoDS_Shape&              after,
                     const std::vector<TopoDS_Shape>& selected,
                     const BlendModified&             modified)
{
  struct Face
  {
    TopoDS_Face          face;
    Bnd_Box              box;
    std::vector<gp_Pnt>  corners;
    bool                 made;
    std::vector<Bnd_Box> patches;
    bool                 tight = false;
  };

  std::vector<Face> faces;
  ShapeSet          supports;
  for (TopExp_Explorer old(before, TopAbs_FACE); old.More(); old.Next())
  {
    supports.Add(old.Current());
    for (const auto& changed : modified(old.Current()))
      if (changed.ShapeType() == TopAbs_FACE)
        supports.Add(changed);
  }
  for (TopExp_Explorer f(after, TopAbs_FACE); f.More(); f.Next())
  {
    Face item{TopoDS::Face(f.Current()), {}, {}, true, {}};
    item.made = !supports.Contains(item.face);
    BRepBndLib::Add(item.face, item.box, false);
    for (TopExp_Explorer v(item.face, TopAbs_VERTEX); v.More(); v.Next())
      item.corners.push_back(BRep_Tool::Pnt(TopoDS::Vertex(v.Current())));
    faces.push_back(std::move(item));
  }
  // Retain remote original faces as obstacles even if the rolling builder
  // accidentally drops their entire inner shell from its result.
  for (TopExp_Explorer old(before, TopAbs_FACE); old.More(); old.Next())
  {
    bool present = false, incident = false;
    for (const auto& f : faces)
      present |= f.face.IsSame(old.Current());
    if (present || !modified(old.Current()).IsEmpty())
      continue;
    Face item{TopoDS::Face(old.Current()), {}, {}, false, {}};
    for (TopExp_Explorer v(old.Current(), TopAbs_VERTEX); v.More(); v.Next())
    {
      const auto p = BRep_Tool::Pnt(TopoDS::Vertex(v.Current()));
      item.corners.push_back(p);
      for (const auto& edge : selected)
        for (TopExp_Explorer ev(edge, TopAbs_VERTEX); ev.More(); ev.Next())
          incident |=
            p.Distance(BRep_Tool::Pnt(TopoDS::Vertex(ev.Current()))) < BRepExt::tolerance::length;
    }
    if (incident)
      continue;
    BRepBndLib::Add(item.face, item.box, false);
    faces.push_back(std::move(item));
  }
  for (size_t i = 0; i < faces.size(); ++i)
    for (size_t j = i + 1; j < faces.size(); ++j)
    {
      auto& a = faces[i];
      auto& b = faces[j];
      if ((!a.made && !b.made) || a.box.IsOut(b.box))
        continue;
      bool adjacent = false;
      for (const auto& p : a.corners)
        for (const auto& q : b.corners)
          adjacent |= p.Distance(q) < BRepExt::tolerance::length;
      if (adjacent)
        continue;
      // The exact box of a face of an elementary kind is dear (some milliseconds for the
      // torus of every rounded rim), and is worked out only for a pair that the loose
      // boxes and the corners have not already told apart.
      const auto tighten = [](Face& item) {
        if (item.tight)
          return;
        item.tight = true;
        if (BRepAdaptor_Surface(item.face).GetType() > GeomAbs_Torus)
          return;
        Bnd_Box exact;
        BRepBndLib::AddOptimal(item.face, exact, false, false);
        item.box = exact;
      };
      tighten(a);
      tighten(b);
      if (a.box.IsOut(b.box))
        continue;
      const auto patches = [](Face& item) {
        if (!item.patches.empty())
          return;
        BRepAdaptor_Surface surface(item.face);
        const auto          type = surface.GetType();
        if (type != GeomAbs_BSplineSurface && type != GeomAbs_SurfaceOfExtrusion)
        {
          item.patches.push_back(item.box);
          return;
        }
        // A closed spline wall's single AABB covers every hole in the cap.
        // Bound its parametric pieces separately; these are enclosing
        // surface boxes, not sampled points or a mesh collision test.
        const double u0 = surface.FirstUParameter(), u1 = surface.LastUParameter();
        const double v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
        if (type == GeomAbs_BSplineSurface)
        {
          for (const auto& patch : GeomAlgoExt::surfacePatches(surface.BSpline(), u0, u1, v0, v1))
            for (int i = 0; i < 4; ++i)
              for (int j = 0; j < 4; ++j)
              {
                auto piece = Handle(Geom_BezierSurface)::DownCast(patch->Copy());
                piece->Segment(i * .25, (i + 1) * .25, j * .25, (j + 1) * .25);
                auto box = GeomAlgoExt::bezierBounds(piece);
                box.Enlarge(BRepExt::tolerance::length);
                item.patches.push_back(box);
              }
        }
        else
          for (int k = 0; k < 32; ++k)
          {
            Bnd_Box box;
            BndLib_AddSurface::Add(surface,
                                   u0 + (u1 - u0) * k / 32,
                                   u0 + (u1 - u0) * (k + 1) / 32,
                                   v0,
                                   v1,
                                   BRepExt::tolerance::length,
                                   box);
            item.patches.push_back(box);
          }
      };
      patches(a);
      patches(b);
      bool overlap = false;
      for (const auto& x : a.patches)
        for (const auto& y : b.patches)
          overlap |= !x.IsOut(y);
      if (!overlap)
        continue;
      BRepExtrema_DistShapeShape distance(a.face, b.face);
      if (distance.IsDone() && distance.Value() <= BRepExt::tolerance::length)
        fail("fillet",
             no_solution,
             "new blend face runs into a nonadjacent face behind its supports");
    }
}
} // namespace BRepExt
