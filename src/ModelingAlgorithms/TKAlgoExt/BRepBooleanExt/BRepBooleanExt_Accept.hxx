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
#include <BRepAlgoExt_Properties.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>

namespace BRepExt
{
// Sheet/curve components are allowed in mixed-dimensional Booleans. Every solid component
// must have closed shells and no free, non-degenerate edges. Empty Boolean results are valid.
inline bool accepts(const TopoDS_Shape& shape, bool solidExpected = false)
{
  if (shape.IsNull() || !BRepCheck_Analyzer(shape).IsValid())
    return false;
  bool solid = false;
  for (TopExp_Explorer body(shape, TopAbs_SOLID); body.More(); body.Next())
  {
    solid = true;
    for (TopExp_Explorer shell(body.Current(), TopAbs_SHELL); shell.More(); shell.Next())
      if (!BRep_Tool::IsClosed(shell.Current()))
        return false;
    // Internal imprint edges can belong to one face without being a free shell edge.
    // Count boundary occurrences, including both uses of a seam, instead of ancestors.
    NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> uses;
    for (TopExp_Explorer face(body.Current(), TopAbs_FACE); face.More(); face.Next())
      for (TopExp_Explorer edge(face.Current(), TopAbs_EDGE); edge.More(); edge.Next())
      {
        const auto& shape = edge.Current();
        if (shape.Orientation() == TopAbs_INTERNAL || shape.Orientation() == TopAbs_EXTERNAL
            || BRep_Tool::Degenerated(TopoDS::Edge(shape)))
          continue;
        if (!uses.IsBound(shape))
          uses.Bind(shape, 0);
        ++uses.ChangeFind(shape);
      }
    for (decltype(uses)::Iterator edge(uses); edge.More(); edge.Next())
      if (edge.Value() < 2)
        return false;
  }
  return !solidExpected || solid;
}

} // namespace BRepExt
