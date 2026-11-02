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

#pragma once
#include <Standard_Macro.hxx>
#include <GeomAlgoExt_Coordinates.hxx>

// Small pieces of geometry and topology that several operations share, each kept in one place.

#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRep_Builder.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>

#include <cmath>
#include <optional>
#include <utility>
#include <vector>

class IntTools_Context;

namespace BRepExt
{

//! Distinct shapes in the order they were added, each found again by its index.
using ShapeSet = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
//! Each shape with the shapes next to it, as TopExp::MapShapesAndAncestors fills it: the faces
//! of each edge, the edges of each vertex.
using Adjacency =
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>;

//! `value` moved by whole periods to the place nearest `toward`.
inline double nearPeriod(double value, double toward, double period)
{
  return value + period * std::round((toward - value) / period);
}

//! The normal of `face` at (u, v) as the face is oriented, outward on a solid; false where the
//! surface has none.
inline bool outwardNormal(const TopoDS_Face& face, double u, double v, gp_Pnt& at, gp_Dir& normal)
{
  const BRepAdaptor_Surface surface(face);
  BRepLProp_SLProps         props(surface, u, v, 1, Precision::Confusion());
  if (!props.IsNormalDefined())
    return false;
  at     = props.Value();
  normal = props.Normal();
  if (face.Orientation() == TopAbs_REVERSED)
    normal.Reverse();
  return true;
}

//! Distances to one shape from many points or shapes. BRepExtrema_DistShapeShape works out the
//! boxes of every vertex, edge and face it measures between when it is first asked; with `to`
//! loaded once as its second shape, those of `to` serve every measure, and each comes out as a
//! fresh BRepExtrema_DistShapeShape(from, to) would.
class DistanceTo
{
public:
  explicit DistanceTo(const TopoDS_Shape& to) { measure_.LoadS2(to); }

  //! The least distance from `from`; none where BRepExtrema finds no solution.
  std::optional<double> from(const TopoDS_Shape& from)
  {
    measure_.LoadS1(from);
    if (!measure_.Perform())
      return std::nullopt;
    return measure_.Value();
  }

  std::optional<double> from(const gp_Pnt& point)
  {
    return from(BRepBuilderAPI_MakeVertex(point).Vertex());
  }

private:
  BRepExtrema_DistShapeShape measure_;
};

//! Parameters inside `face` that a grid of its parameter box may miss altogether (a ring, a
//! slot): the middle of every stretch inside the face of `lines` lines across the box each way,
//! at the middles of its equal columns and rows, as OpenCascade's hatcher finds them. Empty only
//! where the face has no bounded box or the hatcher fails.
Standard_EXPORT std::vector<gp_Pnt2d> parametersAcross(const TopoDS_Face&              face,
                                                       int                             lines,
                                                       const Handle(IntTools_Context)& context);

//! The axis of a planar face's plane, turned round where the face is reversed.
inline gp_Dir orientedPlaneAxis(const TopoDS_Face& face)
{
  gp_Dir axis = BRepAdaptor_Surface(face).Plane().Axis().Direction();
  if (face.Orientation() == TopAbs_REVERSED)
    axis.Reverse();
  return axis;
}

//! The one solid in `shape`, be it wrapped in compounds; a null shape when there is not one.
inline TopoDS_Shape singleSolid(const TopoDS_Shape& shape)
{
  TopExp_Explorer solids(shape, TopAbs_SOLID);
  if (!solids.More())
    return {};
  const TopoDS_Shape solid = solids.Current();
  solids.Next();
  return solids.More() ? TopoDS_Shape() : solid;
}

//! The one solid in `shape` where there is exactly one, else `shape` as it is.
inline TopoDS_Shape solidIfSingle(const TopoDS_Shape& shape)
{
  const TopoDS_Shape solid = singleSolid(shape);
  return solid.IsNull() ? shape : solid;
}

//! An edge that is one point in space, `vertex`, where a surface closes to a pole or an apex: in
//! the surface's parameters it is `pcurve` from `from` to `to`.
inline TopoDS_Edge degenerateEdge(const Handle(Geom2d_Curve)& pcurve,
                                  const Handle(Geom_Surface)& surface,
                                  const TopLoc_Location&      location,
                                  const TopoDS_Vertex&        vertex,
                                  double                      from,
                                  double                      to,
                                  double                      tolerance)
{
  BRep_Builder builder;
  TopoDS_Edge  edge;
  builder.MakeEdge(edge);
  builder.UpdateEdge(edge, pcurve, surface, location, tolerance);
  builder.Degenerated(edge, true);
  builder.Add(edge, vertex.Oriented(TopAbs_FORWARD));
  builder.Add(edge, vertex.Oriented(TopAbs_REVERSED));
  builder.Range(edge, from, to);
  return edge;
}

//! Solve the `n` equations in the first `n` columns of `rows` for the three right-hand sides in
//! the columns after them, by Gauss-Jordan elimination with partial pivoting; the solutions are
//! left in those three columns. False where a pivot is below 1e-14: the system is singular.
template <typename Rows>
bool solveForThree(Rows& rows, int n)
{
  for (int col = 0; col < n; ++col)
  {
    int pivot = col;
    for (int i = col + 1; i < n; ++i)
      if (std::abs(rows[i][col]) > std::abs(rows[pivot][col]))
        pivot = i;
    if (std::abs(rows[pivot][col]) < 1e-14)
      return false;
    std::swap(rows[col], rows[pivot]);
    const double lead = rows[col][col];
    for (int j = col; j < n + 3; ++j)
      rows[col][j] /= lead;
    for (int i = 0; i < n; ++i)
      if (i != col)
      {
        const double f = rows[i][col];
        for (int j = col; j < n + 3; ++j)
          rows[i][j] -= f * rows[col][j];
      }
  }
  return true;
}

} // namespace BRepExt
