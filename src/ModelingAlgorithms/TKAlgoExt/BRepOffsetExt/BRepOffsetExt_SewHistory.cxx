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

#include <BRepOffsetExt_SewHistory.hxx>
#include <BRepAlgoExt_Geometry.hxx>
#include <BRepBooleanExt_Accept.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <cmath>

namespace BRepExt
{
TopoDS_Shape sewn(const BRepBuilderAPI_Sewing& sewing, const TopoDS_Shape& added)
{
  if (sewing.IsModified(added))
    return sewing.Modified(added);
  if (sewing.IsModifiedSubShape(added))
    return sewing.ModifiedSubShape(added);
  return added;
}

std::vector<TopoDS_Shape> sewingHistory(
  const BRepBuilderAPI_Sewing&                            sewing,
  const std::vector<TopoDS_Shape>&                        sources,
  const TopoDS_Shape&                                     result,
  std::initializer_list<TopAbs_ShapeEnum>                 types,
  BRepTools_History&                                      history,
  const std::function<TopoDS_Shape(const TopoDS_Shape&)>& added)
{
  std::vector<TopoDS_Shape> missing;
  for (auto type : types)
  {
    ShapeSet delivered;
    TopExp::MapShapes(result, type, delivered);
    for (const auto& owner : sources)
    {
      ShapeSet entities;
      TopExp::MapShapes(owner, type, entities);
      for (const auto& source : entities)
      {
        const auto successor = sewn(sewing, added ? added(source) : source);
        bool       found     = false;
        if (!successor.IsNull())
        {
          ShapeSet products;
          TopExp::MapShapes(successor, type, products);
          for (const auto& product : products)
            if (delivered.Contains(product))
            {
              history.AddModified(source, product);
              found = true;
            }
        }
        if (!found)
          missing.push_back(source);
      }
    }
  }
  return missing;
}

TopoDS_Shape closedSolid(const char*                  op,
                         const BRepBuilderAPI_Sewing& sewing,
                         const char*                  reason,
                         double                       minimumVolume)
{
  const auto& shell = sewing.SewedShape();
  if (shell.IsNull() || shell.ShapeType() != TopAbs_SHELL || sewing.NbFreeEdges() != 0)
    return {};
  auto solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shell)).Solid();
  if (!BRepLib::OrientClosedSolid(solid))
    fail(op, reason, "The closed shell could not be oriented.");
  if (!(estimatedVolume(solid) > minimumVolume))
    fail(op, reason, "The closed shell has no positive volume.");
  return solid;
}
} // namespace BRepExt
