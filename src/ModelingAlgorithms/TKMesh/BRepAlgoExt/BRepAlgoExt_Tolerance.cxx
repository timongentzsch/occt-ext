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

#include <BRepAlgoExt_Tolerance.hxx>

#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <cmath>
#include <limits>

namespace BRepExt::tolerance
{

double own(const TopoDS_Shape& shape)
{
  switch (shape.ShapeType())
  {
    case TopAbs_VERTEX:
      return BRep_Tool::Tolerance(TopoDS::Vertex(shape));
    case TopAbs_EDGE:
      return BRep_Tool::Tolerance(TopoDS::Edge(shape));
    case TopAbs_FACE:
      return BRep_Tool::Tolerance(TopoDS::Face(shape));
    default:
      return 0;
  }
}

double largest(const TopoDS_Shape& shape)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> seen;
  TopExp::MapShapes(shape, seen);
  double value = 0;
  for (int i = 1; i <= seen.Extent(); ++i)
  {
    const double local = own(seen(i));
    if (!std::isfinite(local))
      return std::numeric_limits<double>::infinity();
    value = std::max(value, local);
  }
  return value;
}

} // namespace BRepExt::tolerance
