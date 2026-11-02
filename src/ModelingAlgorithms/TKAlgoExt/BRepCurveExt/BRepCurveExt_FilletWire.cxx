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

// Rounding the corners of a planar wire.
//
// A circle of radius r tangent to the two pieces that meet at a corner has its centre on the
// offset of each piece by r, on the side the wire turns to; the tangent points are the feet of the
// centre on the pieces. For two lines this is the closed form w = r / tan(theta / 2) of the
// distance from the corner to each tangent point (theta the interior angle); for curved pieces it
// is the root of a one-dimensional distance function that OCCT's planar fillet algorithm finds by
// Newton's method from several starts (ChFi2d_FilletAPI), so the pieces are trimmed exactly at the
// tangent points and the arc is exact. A negative radius cuts the corner along the chord between
// the same two tangent points. Corners are made in the order the wire runs; a piece that two
// selected corners share is trimmed at both ends in turn.

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepLib_FindSurface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <ChFi2d_FilletAPI.hxx>
#include <Geom_Plane.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pln.hxx>

#include <cmath>
#include <string>
#include <vector>

namespace BRepExt
{

EvolutionData filletWire(OcctKernel&                  kernel,
                         uint32_t                     wireId,
                         const std::vector<uint32_t>& vertexIds,
                         double                       radius,
                         double                       tolerance)
{
  const char* const   op    = "filletWire";
  const TopoDS_Shape& shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE)
    fail(op, invalid_input, "a wire is required");
  if (!std::isfinite(radius) || radius == 0.0 || !std::isfinite(tolerance) || tolerance <= 0.0)
    fail(op, invalid_input, "a nonzero finite radius and a positive tolerance are required");
  if (vertexIds.empty())
    fail(op, invalid_input, "no vertex was chosen");
  const TopoDS_Wire wire = TopoDS::Wire(shape);

  const std::vector<TopoDS_Shape> chosen = subShapes(kernel, op, wire, vertexIds, TopAbs_VERTEX);
  for (size_t i = 0; i < chosen.size(); ++i)
    for (size_t k = i + 1; k < chosen.size(); ++k)
      if (chosen[i].IsSame(chosen[k]))
        fail(op, invalid_input, "a vertex is given twice");

  // The plane every piece lies in: the circles are in it, and the algorithm that finds them works
  // in it.
  const BRepLib_FindSurface plane(wire, tolerance, true);
  if (!plane.Found())
    fail(op, unsupported_geometry, "the wire is not planar");
  const gp_Pln pln = Handle(Geom_Plane)::DownCast(plane.Surface())->Pln();

  // The pieces in the order the wire runs, and the vertex that follows each.
  std::vector<TopoDS_Edge>   pieces;
  std::vector<TopoDS_Vertex> after;
  for (BRepTools_WireExplorer it(wire); it.More(); it.Next())
  {
    pieces.push_back(it.Current());
    after.push_back(it.CurrentVertex());
  }
  if (pieces.empty())
    fail(op, invalid_input, "the wire has no edges");
  const size_t n      = pieces.size();
  const bool   closed = BRep_Tool::IsClosed(wire);
  // The wire explorer gives the vertex a piece starts at; the one it ends at starts the next.
  const auto vertexAfter = [&](size_t i) { return after[(i + 1) % n]; };

  std::vector<std::vector<TopoDS_Edge>> beforeNext(n); // what is inserted after piece i
  std::vector<bool>                     used(chosen.size(), false);
  for (size_t i = 0; i < n; ++i)
  {
    if (!closed && i + 1 == n)
      break; // the end of an open wire is not a corner
    const TopoDS_Vertex corner = vertexAfter(i);
    int                 which  = -1;
    for (size_t k = 0; k < chosen.size(); ++k)
      if (chosen[k].IsSame(corner))
        which = int(k);
    if (which < 0)
      continue;
    used[size_t(which)]   = true;
    const size_t     next = (i + 1) % n;
    ChFi2d_FilletAPI maker(pieces[i], pieces[next], pln);
    maker.Perform(std::abs(radius));
    const gp_Pnt at = BRep_Tool::Pnt(corner);
    TopoDS_Edge  arc, first, second;
    if (maker.NbResults(at) < 1)
      fail(op,
           no_solution,
           "a corner is too tight for the radius: the circle does not fit between its pieces");
    arc = maker.Result(at, first, second);
    if (arc.IsNull() || first.IsNull() || second.IsNull())
      fail(op, no_solution, "a corner could not be rounded");
    pieces[i]    = first;
    pieces[next] = second;
    if (radius < 0.0)
    {
      // The chord between the arc's ends, the same tangent points.
      TopoDS_Vertex a, b;
      TopExp::Vertices(arc, a, b, true);
      arc = BRepBuilderAPI_MakeEdge(BRep_Tool::Pnt(a), BRep_Tool::Pnt(b)).Edge();
    }
    beforeNext[i].push_back(arc);
  }
  for (size_t k = 0; k < chosen.size(); ++k)
    if (!used[k])
      fail(op, invalid_input, "a vertex is not a corner of the wire");

  BRepBuilderAPI_MakeWire made;
  for (size_t i = 0; i < n; ++i)
  {
    made.Add(pieces[i]);
    for (const TopoDS_Edge& arc : beforeNext[i])
      made.Add(arc);
  }
  if (!made.IsDone())
    fail(op, no_solution, "the trimmed pieces and the arcs do not join up");
  const TopoDS_Wire result = made.Wire();
  const std::string reason = invalidity(result);
  if (!reason.empty())
    fail(op, no_solution, "the rounded wire is not valid");
  return stored(kernel, result);
}

} // namespace BRepExt
