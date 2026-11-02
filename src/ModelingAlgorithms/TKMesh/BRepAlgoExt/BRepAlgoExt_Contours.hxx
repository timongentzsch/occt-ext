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

#include <BRepAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Precision.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <vector>

namespace BRepExt
{

// Read an already seeded maker's own tangent contours without building a body.
// Keep the explicit selection first, including any edge the maker did not adopt.
// Both native construction and medial preflight use this same propagation census.
template <typename Maker>
std::vector<TopoDS_Shape> blendContourMembers(Maker&                           maker,
                                              const std::vector<TopoDS_Shape>& selected)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> runs;
  for (const auto& edge : selected)
    runs.Add(edge);
  for (int contour = 1; contour <= maker.NbContours(); ++contour)
    for (int i = 1; i <= maker.NbEdges(contour); ++i)
      runs.Add(maker.Edge(contour, i));
  std::vector<TopoDS_Shape> result;
  for (int i = 1; i <= runs.Extent(); ++i)
    result.push_back(runs(i));
  return result;
}

// GeomAPI_Interpolate's periodic spline is only C1 at its closing knot. OCCT's
// closed-contour walker can miss that knot when seeded on another trimmed part
// of the same curve: merely reversing the selected edges changes the blend.
// Seed such a contour at the basis curve's closing knot. This changes neither
// geometry nor topology, and also handles selecting only the other half.
// Maker is MakeFillet or MakeChamfer; discovery uses their own tangent chaining.
template <typename Maker>
std::vector<TopoDS_Shape> blendContours(const TopoDS_Shape&              shape,
                                        const std::vector<TopoDS_Shape>& edges)
{
  bool periodic = false;
  for (const auto& edge : edges)
  {
    BRepAdaptor_Curve curve(TopoDS::Edge(edge));
    periodic |= curve.GetType() == GeomAbs_BSplineCurve && curve.IsPeriodic();
  }
  if (!periodic)
    return edges;

  Maker                     contours(shape);
  std::vector<TopoDS_Shape> seeds;
  for (const auto& selected : edges)
  {
    const auto edge = TopoDS::Edge(selected);
    if (contours.Contour(edge))
      continue;
    try
    {
      contours.Add(edge);
    }
    catch (const Standard_Failure&)
    {
      return edges;
    }
    const int    contour = contours.Contour(edge);
    TopoDS_Shape seed    = edge;
    if (contour && contours.ClosedAndTangent(contour))
    {
      for (int i = 1; i <= contours.NbEdges(contour); ++i)
      {
        const auto        candidate = contours.Edge(contour, i);
        BRepAdaptor_Curve curve(candidate);
        if (curve.GetType() != GeomAbs_BSplineCurve || !curve.IsPeriodic())
          continue;
        const auto basis = curve.BSpline();
        if (std::abs(curve.FirstParameter() - basis->FirstParameter()) <= Precision::PConfusion())
        {
          seed = candidate;
          break;
        }
      }
    }
    seeds.push_back(seed);
  }
  return seeds;
}

} // namespace BRepExt
