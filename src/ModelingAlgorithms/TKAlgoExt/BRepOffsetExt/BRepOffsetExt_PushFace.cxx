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

// Push-pull through changes of topology, after Zou and Feng (Push-pull direct modeling of
// solid CAD models, 2019): the solid that the moved faces sweep is built by itself and joined
// to the body (a pull) or taken from it (a push) by a Boolean. Whatever the move does to the
// rest of the body (a step that goes, a hole that breaks out, a floor that goes through the
// far side) is then the Boolean's to resolve, and no case of its own.
//
// The sides of the swept solid have to be faces of the body, extended, or the body gets
// steps. Two swept solids are built here whose sides are exact:
//
//   - flat faces that move one way: the prism of the faces and of the blends that are
//     tangent to them, along the move. A blend moves with its face, as a blend does; a
//     neighbour that the move runs along is extended by the prism's side.
//   - any other face: the solid between the face and its offset, with sides along the
//     normals. Where a neighbour is tangent that side is a step face, which is what is asked
//     for there; a sharp neighbour has to run along the normals.
//
// A sharp neighbour that slants is neither: it is refused, naming the face.

#include <BRepOffsetExt_Reconstruct.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepOffsetExt_Solid.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffset_MakeOffset.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>

#include <BRepClass_FaceClassifier.hxx>
#include <cmath>
#include <deque>

namespace BRepExt
{

namespace
{

//! The outward normals of `face` at five places along `edge`.
std::vector<gp_Dir> normalsAlong(const TopoDS_Edge& edge, const TopoDS_Face& face)
{
  std::vector<gp_Dir> out;
  BRepAdaptor_Curve2d onFace(edge, face);
  for (int i = 0; i < 5; i++)
  {
    const double t =
      onFace.FirstParameter() + (onFace.LastParameter() - onFace.FirstParameter()) * (i + 0.5) / 5;
    const gp_Pnt2d uv = onFace.Value(t);
    gp_Pnt         at;
    gp_Dir         normal;
    if (outwardNormal(face, uv.X(), uv.Y(), at, normal))
      out.push_back(normal);
  }
  return out;
}

//! The least and the greatest of the face's outward normal along `direction`.
std::pair<double, double> leanOf(const TopoDS_Face& face, const gp_Dir& direction)
{
  const BRepAdaptor_Surface surface(face);
  double                    low = 1.0, high = -1.0;
  for (int i = 0; i <= 8; i++)
    for (int j = 0; j <= 8; j++)
    {
      gp_Pnt       at;
      gp_Dir       normal;
      const double u =
        surface.FirstUParameter() + (surface.LastUParameter() - surface.FirstUParameter()) * i / 8;
      const double v =
        surface.FirstVParameter() + (surface.LastVParameter() - surface.FirstVParameter()) * j / 8;
      // A sample off the trimmed face (a triangular patch fills a part of its chart) is
      // not the face.
      if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-9).State() == TopAbs_OUT
          || !outwardNormal(face, u, v, at, normal))
        continue;
      low  = std::min(low, normal.Dot(direction));
      high = std::max(high, normal.Dot(direction));
    }
  return {low, high};
}

} // namespace

TopoDS_Shape sweepFaces(const TopoDS_Shape&              solid,
                        const std::vector<TopoDS_Shape>& moved,
                        double                           distance,
                        double                           tolerance,
                        Lineage&                         faces,
                        std::string&                     why,
                        std::vector<TopoDS_Shape>&       at)
{
  Adjacency facesOf;
  TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, facesOf);

  // Flat faces that all face one way move as a prism.
  bool   flat = true;
  gp_Dir direction;
  for (size_t i = 0; i < moved.size() && flat; i++)
  {
    const TopoDS_Face         face = TopoDS::Face(moved[i]);
    const BRepAdaptor_Surface surface(face);
    gp_Pnt                    point;
    gp_Dir                    normal;
    flat =
      surface.GetType() == GeomAbs_Plane
      && outwardNormal(face, surface.FirstUParameter(), surface.FirstVParameter(), point, normal)
      && (i == 0 || normal.Dot(direction) > 1.0 - alongSlack);
    if (i == 0)
      direction = normal;
  }

  // The faces that move: with flat faces, the blends tangent to them as well.
  FaceSet                   chain;
  std::vector<TopoDS_Shape> order;
  std::deque<TopoDS_Shape>  open;
  for (const TopoDS_Shape& face : moved)
    if (chain.Add(face))
    {
      order.push_back(face);
      open.push_back(face);
    }
  while (flat && !open.empty())
  {
    const TopoDS_Face face = TopoDS::Face(open.front());
    open.pop_front();
    for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
      for (const TopoDS_Shape& other : facesOf.FindFromKey(edges.Current()))
      {
        if (chain.Contains(other)
            || !tangentAlong(TopoDS::Edge(edges.Current()), face, TopoDS::Face(other)))
          continue;
        const auto [low, high] = leanOf(TopoDS::Face(other), direction);
        // The move runs along it: it stays, and the prism's side extends it.
        if (std::max(std::abs(low), std::abs(high)) <= alongSlack)
          continue;
        if (low < -alongSlack)
        {
          why = "a blend beside the moved face turns away from the move further on; it cannot move "
                "with the face";
          at  = {other};
          return {};
        }
        chain.Add(other);
        order.push_back(other);
        open.push_back(other);
      }
  }

  // The sides of the swept solid: each must be a face of the body, or a step at a tangent one.
  for (const TopoDS_Shape& face : order)
    for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
    {
      const TopoDS_Edge edge = TopoDS::Edge(edges.Current());
      if (BRep_Tool::Degenerated(edge))
        continue;
      for (const TopoDS_Shape& other : facesOf.FindFromKey(edge))
      {
        if (chain.Contains(other) || tangentAlong(edge, TopoDS::Face(face), TopoDS::Face(other)))
          continue;
        const std::vector<gp_Dir> theirs = normalsAlong(edge, TopoDS::Face(other));
        const std::vector<gp_Dir> ours   = normalsAlong(edge, TopoDS::Face(face));
        for (size_t i = 0; i < theirs.size() && i < ours.size(); i++)
          if (std::abs(theirs[i].Dot(flat ? direction : ours[i])) > alongSlack)
          {
            why = "a face beside the moved one slants against the move, and the move changes "
                  "which faces the body has; such a neighbour cannot be extended through that";
            at  = {other};
            return {};
          }
      }
    }

  // The swept solid, and where in it each moved face ends up. It is swept from a copy of the
  // faces: a prism shares the edges of what it is made from and adds its sides' curves to
  // them, and the body's own edges are not to change.
  TopoDS_Shell original;
  BRep_Builder().MakeShell(original);
  for (const TopoDS_Shape& face : order)
    BRep_Builder().Add(original, face);
  BRepBuilderAPI_Copy       copy(original);
  const TopoDS_Shape        skin = copy.Shape();
  std::vector<TopoDS_Shape> copies;
  for (const TopoDS_Shape& face : order)
    copies.push_back(copy.ModifiedShape(face));
  TopoDS_Shape                           tool;
  std::vector<std::vector<TopoDS_Shape>> ends(order.size());
  try
  {
    if (flat)
    {
      BRepPrimAPI_MakePrism prism(skin, gp_Vec(direction) * distance);
      if (!prism.IsDone())
        return {};
      tool = prism.Shape();
      for (size_t i = 0; i < order.size(); i++)
        ends[i] = {prism.LastShape(copies[i])};
    }
    else
    {
      // A face that the distance turns inside out has no offset to sweep to.
      const OffsetFate fate = classifyForOffset(skin, distance);
      if (!fate.vanish.empty() || !fate.partial.empty())
      {
        why = "the distance exceeds the radius of curvature of a moved face";
        at  = fate.vanish.empty() ? fate.partial : fate.vanish;
        return {};
      }
      BRepOffset_MakeOffset thick;
      thick.Initialize(skin,
                       distance,
                       tolerance,
                       BRepOffset_Skin,
                       false,
                       false,
                       GeomAbs_Intersection,
                       true);
      thick.MakeOffsetShape();
      if (!thick.IsDone() || thick.Shape().IsNull())
        return {};
      tool = thick.Shape();
      for (size_t i = 0; i < order.size(); i++)
        for (const TopoDS_Shape& image : thick.Generated(copies[i]))
          if (image.ShapeType() == TopAbs_FACE)
            ends[i].push_back(image);
    }
  }
  catch (const Standard_Failure&)
  {
    return {};
  }
  NCollection_List<TopoDS_Shape> bodies, tools;
  bodies.Append(solid);
  for (TopExp_Explorer solids(tool, TopAbs_SOLID); solids.More(); solids.Next())
    tools.Append(solids.Current());
  if (tools.IsEmpty() || !BRepCheck_Analyzer(tool).IsValid())
    return {};

  BRepAlgoAPI_Fuse              join;
  BRepAlgoAPI_Cut               take;
  BRepAlgoAPI_BooleanOperation& boolean = distance > 0
                                            ? static_cast<BRepAlgoAPI_BooleanOperation&>(join)
                                            : static_cast<BRepAlgoAPI_BooleanOperation&>(take);
  boolean.SetArguments(bodies);
  boolean.SetTools(tools);
  boolean.SetNonDestructive(true);
  boolean.Build();
  if (!boolean.IsDone() || boolean.HasErrors() || boolean.Shape().IsNull())
    return {};

  // The prism's sides lie in the neighbours they extend: one face each again. Edges that
  // the body had between faces that stay are kept, whatever lies beside them.
  ShapeUpgrade_UnifySameDomain                           unify(boolean.Shape(), true, true, false);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> keep;
  for (int i = 1; i <= facesOf.Extent(); i++)
  {
    bool ofChain = false;
    for (const TopoDS_Shape& face : facesOf(i))
      ofChain = ofChain || chain.Contains(face);
    if (!ofChain)
      keep.Add(facesOf.FindKey(i));
  }
  unify.KeepShapes(keep);
  unify.Build();
  TopoDS_Shape result = unify.Shape();
  int          count  = 0;
  for (TopExp_Explorer solids(result, TopAbs_SOLID); solids.More(); solids.Next())
    count++;
  if (count == 0 || !BRepCheck_Analyzer(result).IsValid())
  {
    why = "nothing of the body is left, or what is left is not a valid solid";
    return {};
  }
  if (count == 1)
    result = singleSolid(result);
  // A pull adds and a push takes away.
  if ((estimatedVolume(result) - estimatedVolume(solid)) * distance <= 0.0)
  {
    why = "the move changed nothing";
    return {};
  }

  faces.through(boolean);
  faces.throughHistory(unify.History());
  Lineage swept(tool);
  swept.through(boolean);
  swept.throughHistory(unify.History());
  swept.within(result);
  for (size_t i = 0; i < order.size(); i++)
  {
    std::vector<TopoDS_Shape> now;
    for (const TopoDS_Shape& end : ends[i])
      for (const TopoDS_Shape& successor : swept.of(end))
        now.push_back(successor);
    // Apart from the Boolean's word: where the face is now, it is the distance from where
    // it was. A flat one and its blends are moved back and must lie on the old surface; any
    // other must be the distance from it.
    const Handle(Geom_Surface) old = BRep_Tool::Surface(TopoDS::Face(order[i]));
    for (const TopoDS_Shape& face : now)
      for (gp_Pnt point : pointsInside(TopoDS::Face(face), 3))
      {
        if (flat)
          point.Translate(gp_Vec(direction) * -distance);
        const GeomAPI_ProjectPointOnSurf onOld(point, old);
        if (onOld.NbPoints() == 0
            || std::abs(onOld.LowerDistance() - (flat ? 0.0 : std::abs(distance)))
                 > clearanceSlack * tolerance)
        {
          why = "a moved face did not come to lie at the distance";
          at  = {order[i]};
          return {};
        }
      }
    faces.set(order[i], now);
  }
  // A face that went through the far side of the body is gone, and the body is open there.
  faces.within(result);
  return result;
}

} // namespace BRepExt
