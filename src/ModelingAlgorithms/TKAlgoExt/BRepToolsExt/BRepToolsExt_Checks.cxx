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

#include <BRepToolsExt_Checks.hxx>

#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Shell.hxx>
#include <BRepCheck_Solid.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <stdexcept>
#include <unordered_set>

namespace BRepExt
{

std::vector<TopoDS_Shape> subShapes(const OcctKernel&            kernel,
                                    const char*                  operation,
                                    const TopoDS_Shape&          owner,
                                    const std::vector<uint32_t>& ids,
                                    TopAbs_ShapeEnum             type)
{
  const char* noun = type == TopAbs_FACE ? "face" : type == TopAbs_VERTEX ? "vertex" : "edge";
  if (ids.empty())
    fail(operation, invalid_input, std::string("no ") + noun + " given");

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> owned;
  TopExp::MapShapes(owner, type, owned);

  std::vector<TopoDS_Shape> found;
  found.reserve(ids.size());
  std::unordered_set<int> seen;
  for (uint32_t id : ids)
  {
    const int index = owned.FindIndex(get(kernel, id));
    if (index == 0)
    {
      fail(operation,
           invalid_input,
           std::string("selection is not ") + (type == TopAbs_FACE ? "a " : "an ") + noun
             + " of the shape being modified");
    }
    if (!seen.insert(index).second)
    {
      fail(operation, invalid_input, std::string(noun) + " is given twice");
    }
    found.push_back(owned.FindKey(index));
  }
  return found;
}

// The statuses that an analyzer of `shape` holds, as "kind:status" words.
static std::string describe(const BRepCheck_Analyzer& check, const TopoDS_Shape& shape)
{
  std::string what;
  for (const auto kind :
       {TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SHELL, TopAbs_SOLID})
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> all;
    TopExp::MapShapes(shape, kind, all);
    for (int k = 1; k <= all.Extent(); ++k)
    {
      const auto& found = check.Result(all(k));
      if (found.IsNull())
        continue;
      for (const auto status : found->Status())
        if (status != BRepCheck_NoError && what.size() < 300)
          what += " " + std::to_string(int(kind)) + ":" + std::to_string(int(status));
      for (found->InitContextIterator(); found->MoreShapeInContext(); found->NextShapeInContext())
        for (const auto status : found->StatusOnShape())
          if (status != BRepCheck_NoError && what.size() < 300)
            what += " " + std::to_string(int(kind)) + "in"
                    + std::to_string(int(found->ContextualShape().ShapeType())) + ":"
                    + std::to_string(int(status));
    }
  }
  return what;
}

std::string invalidity(const TopoDS_Shape& shape)
{
  const BRepCheck_Analyzer check(shape);
  if (check.IsValid())
    return {};
  const std::string what = describe(check, shape);
  return what.empty() ? " unknown" : what;
}

std::string invalidity(const TopoDS_Shape& shape, const TopoDS_Shape& earlier)
{
  if (TopExp_Explorer(shape, TopAbs_SOLID).More() == false)
    return invalidity(shape);
  // A face that the operation did not touch is the earlier body's own face, with its wires, edges
  // and vertices as they were: it keeps the verdict it had, and is not checked again. The faces
  // that are new are checked with everything under them, then each solid and shell as a whole
  // (connected, closed, consistently oriented), which is what remains of the full check.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> known, faces;
  TopExp::MapShapes(earlier, TopAbs_FACE, known);
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  std::string what;
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    if (known.Contains(faces(i)))
      continue;
    const BRepCheck_Analyzer check(faces(i));
    if (check.IsValid())
      continue;
    const std::string found = describe(check, faces(i));
    what += found.empty() ? " unknown" : found;
  }
  for (TopExp_Explorer solids(shape, TopAbs_SOLID); solids.More(); solids.Next())
  {
    const TopoDS_Solid solid = TopoDS::Solid(solids.Current());
    BRepCheck_Solid    whole(solid);
    whole.Minimum();
    for (const auto status : whole.Status())
      if (status != BRepCheck_NoError)
        what += " " + std::to_string(int(TopAbs_SOLID)) + ":" + std::to_string(int(status));
    for (TopoDS_Iterator shells(solid); shells.More(); shells.Next())
    {
      if (shells.Value().ShapeType() != TopAbs_SHELL)
        continue;
      BRepCheck_Shell shell(TopoDS::Shell(shells.Value()));
      shell.Minimum();
      shell.InContext(solid);
      for (const auto status : shell.Status())
        if (status != BRepCheck_NoError)
          what += " " + std::to_string(int(TopAbs_SHELL)) + ":" + std::to_string(int(status));
      for (const auto status : shell.StatusOnShape(solid))
        if (status != BRepCheck_NoError)
          what += " " + std::to_string(int(TopAbs_SHELL)) + "in" + std::to_string(int(TopAbs_SOLID))
                  + ":" + std::to_string(int(status));
    }
  }
  return what;
}

} // namespace BRepExt
