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

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Circ.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <GeomLib.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BoundedSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <cmath>
#include <gp_Pln.hxx>

namespace BRepExt
{
namespace
{
TopoDS_Face faceInput(OcctKernel& kernel, const char* op, uint32_t id, double tolerance)
{
  const auto& shape = get(kernel, id);
  if (shape.ShapeType() != TopAbs_FACE || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "a face and positive finite tolerance are required");
  return TopoDS::Face(shape);
}

EvolutionData checked(OcctKernel& kernel, const char* op, const TopoDS_Shape& shape)
{
  if (shape.IsNull() || !BRepCheck_Analyzer(shape).IsValid())
    fail(op, no_solution, "invalid sheet result");
  return stored(kernel, shape);
}
} // namespace

EvolutionData untrimFace(OcctKernel& kernel, uint32_t faceId, bool keepEdges, double tolerance)
{
  const char* op      = "untrimFace";
  const auto  face    = faceInput(kernel, op, faceId, tolerance);
  auto        surface = BRep_Tool::Surface(face);
  while (auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface))
    surface = trimmed->BasisSurface();
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  double a, b, c, d;
  BRepTools::UVBounds(face, a, b, c, d);
  // Infinite analytic supports have no recoverable original extent. Keep their
  // finite trimmed UV envelope; finite NURBS recover their entire support
  // domain.
  if (Precision::IsInfinite(u0) || Precision::IsInfinite(u1))
  {
    u0 = a;
    u1 = b;
  }
  if (Precision::IsInfinite(v0) || Precision::IsInfinite(v1))
  {
    v0 = c;
    v1 = d;
  }
  auto result = BRepBuilderAPI_MakeFace(surface, u0, u1, v0, v1, tolerance).Face();
  result.Orientation(face.Orientation());
  if (!keepEdges)
    return checked(kernel, op, result);
  BRepAlgoAPI_Splitter           split;
  NCollection_List<TopoDS_Shape> arguments, tools;
  arguments.Append(result);
  for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    tools.Append(it.Current());
  split.SetArguments(arguments);
  split.SetTools(tools);
  split.SetNonDestructive(true);
  split.Build();
  if (!split.IsDone())
    fail(op, no_solution, "could not retain original trim edges");
  return checked(kernel, op, split.Shape());
}

EvolutionData isoparam(OcctKernel&                kernel,
                       uint32_t                   faceId,
                       bool                       alongU,
                       const std::vector<double>& parameters,
                       double                     tolerance)
{
  const char* op      = "isoparam";
  const auto  face    = faceInput(kernel, op, faceId, tolerance);
  const auto  surface = BRep_Tool::Surface(face);
  double      u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  BRep_Builder    builder;
  TopoDS_Compound result;
  builder.MakeCompound(result);
  int count = 0;
  for (double f : parameters)
  {
    if (!std::isfinite(f) || f < 0 || f > 1)
      fail(op, invalid_input, "parameters must be normalized to [0,1]");
    auto curve = alongU ? surface->UIso(u0 + f * (u1 - u0)) : surface->VIso(v0 + f * (v1 - v0));
    auto edge  = BRepBuilderAPI_MakeEdge(curve, alongU ? v0 : u0, alongU ? v1 : u1).Edge();
    Untouched<BRepAlgoAPI_Common> clip(edge, face);
    if (!clip.IsDone())
      fail(op, no_solution, "isoparam trim failed");
    for (TopExp_Explorer it(clip.Shape(), TopAbs_EDGE); it.More(); it.Next())
    {
      builder.Add(result, it.Current());
      ++count;
    }
  }
  if (!count)
    fail(op, no_solution, "no isoparam intersects the trimmed face");
  return checked(kernel, op, result);
}

EvolutionData extendSheet(OcctKernel& kernel,
                          uint32_t    faceId,
                          int         side,
                          double      distance,
                          int         continuity,
                          double      tolerance)
{
  const char* op   = "extendSheet";
  const auto  face = faceInput(kernel, op, faceId, tolerance);
  if (side < 0 || side > 3 || !std::isfinite(distance) || distance <= tolerance || continuity < 1
      || continuity > 2)
    fail(op, invalid_input, "a UV boundary, positive distance and continuity 1 or 2 are required");
  int wires = 0;
  for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
    ++wires;
  if (wires != 1)
    fail(op, unsupported_geometry, "extend requires an untrimmed rectangular sheet without holes");
  auto   surface = BRep_Tool::Surface(face);
  double u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  // Verify the face really is its UV rectangle, rather than quietly discarding
  // trims.
  const auto rectangle = BRepBuilderAPI_MakeFace(surface, u0, u1, v0, v1, tolerance).Face();
  Untouched<BRepAlgoAPI_Common> common(face, rectangle);
  if (!common.IsDone())
    fail(op, no_solution, "could not verify rectangular support");
  const auto rectangleId = store(kernel, rectangle);
  if (std::abs(surfaceArea(kernel, faceId) - surfaceArea(kernel, rectangleId))
      > tolerance * std::max(1.0, surfaceArea(kernel, faceId)))
    fail(op, unsupported_geometry, "extend requires an untrimmed rectangular sheet");
  Handle(Geom_BoundedSurface) bounded = GeomConvert::SurfaceToBSplineSurface(
    new Geom_RectangularTrimmedSurface(surface, u0, u1, v0, v1));
  GeomLib::ExtendSurfByLength(bounded, distance, continuity, side < 2, side % 2 == 1);
  auto result = BRepBuilderAPI_MakeFace(bounded, tolerance).Face();
  result.Orientation(face.Orientation());
  return checked(kernel, op, result);
}

std::vector<int> boundaryEdgeIndices(OcctKernel& kernel, uint32_t id)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    adjacency;
  TopExp::MapShapesAndAncestors(get(kernel, id), TopAbs_EDGE, TopAbs_FACE, adjacency);
  std::vector<int> out;
  for (int i = 1; i <= adjacency.Extent(); ++i)
  {
    const auto  e     = TopoDS::Edge(adjacency.FindKey(i));
    const auto& faces = adjacency.FindFromIndex(i);
    if (faces.Extent() == 1 && !BRep_Tool::Degenerated(e)
        && !BRep_Tool::IsClosed(e, TopoDS::Face(faces.First())))
      out.push_back(i - 1);
  }
  return out;
}

EvolutionData offsetSheetEdge(OcctKernel&                  kernel,
                              uint32_t                     faceId,
                              const std::vector<uint32_t>& edgeIds,
                              double                       distance,
                              double                       tolerance)
{
  const char* op   = "offsetSheetEdge";
  const auto  face = faceInput(kernel, op, faceId, tolerance);
  if (!std::isfinite(distance) || std::abs(distance) <= tolerance)
    fail(op, invalid_input, "offset must be nonzero");
  BRepAdaptor_Surface surface(face);
  if (surface.GetType() == GeomAbs_Cylinder)
  {
    const auto      cylinder = surface.Cylinder();
    BRep_Builder    builder;
    TopoDS_Compound result;
    builder.MakeCompound(result);
    const double low = surface.FirstVParameter(), high = surface.LastVParameter();
    const auto   selected = subShapes(kernel, op, face, edgeIds, TopAbs_EDGE);
    if (high - low + selected.size() * distance <= tolerance)
      fail(op, no_solution, "offset collapses the cylindrical band");
    for (const auto& shape : selected)
    {
      const auto        edge = TopoDS::Edge(shape);
      BRepAdaptor_Curve c(edge);
      if (c.GetType() != GeomAbs_Circle
          || std::abs(c.LastParameter() - c.FirstParameter() - 2 * M_PI)
               > BRepExt::tolerance::angularResolution
          || std::abs(c.Circle().Axis().Direction().Dot(cylinder.Axis().Direction()))
               < 1 - BRepExt::tolerance::angle)
        fail(op,
             unsupported_geometry,
             "cylindrical geodesic offset supports complete circular boundary loops only");
      const double v =
        gp_Vec(cylinder.Location(), c.Circle().Location()).Dot(gp_Vec(cylinder.Axis().Direction()));
      double sign;
      if (std::abs(v - low) < tolerance)
        sign = -1;
      else if (std::abs(v - high) < tolerance)
        sign = 1;
      else
        fail(op, unsupported_geometry, "circle must bound the cylindrical band");
      gp_Trsf move;
      move.SetTranslation(gp_Vec(cylinder.Axis().Direction()) * (sign * distance));
      auto moved = BRepBuilderAPI_Transform(edge, move, true).Shape();
      builder.Add(result, BRepBuilderAPI_MakeWire(TopoDS::Edge(moved)).Wire());
    }
    return checked(kernel, op, result);
  }
  if (surface.GetType() != GeomAbs_Plane)
    fail(op,
         unsupported_geometry,
         "geodesic offset currently supports planes and cylindrical circle boundaries");
  BRepBuilderAPI_MakeWire wire;
  for (const auto& edge : subShapes(kernel, op, face, edgeIds, TopAbs_EDGE))
    wire.Add(TopoDS::Edge(edge));
  if (!wire.IsDone())
    fail(op, invalid_input, "selected edges must form one connected wire");
  BRepBuilderAPI_MakeFace  support(surface.Plane(), wire.Wire(), false);
  BRepOffsetAPI_MakeOffset offset(support.Face(), GeomAbs_Intersection, true);
  offset.Perform(distance);
  if (!offset.IsDone())
    fail(op, no_solution, "edge offset failed");
  return checked(kernel, op, offset.Shape());
}
} // namespace BRepExt
