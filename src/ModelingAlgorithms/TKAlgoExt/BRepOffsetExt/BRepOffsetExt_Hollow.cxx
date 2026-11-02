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

// A closed hollow: the body minus its inside at the wall's distance (BRepOffsetExt_Solid.cxx), the
// outer skin untouched.

#include <BRepOffsetExt_Solid.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepFilletExt_Untouched.hxx>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <GeomAbs_JoinType.hxx>

#include <cmath>

namespace BRepExt
{

EvolutionData hollow(OcctKernel&             kernel,
                     uint32_t                shapeId,
                     double                  thickness,
                     double                  tolerance,
                     const std::vector<int>& inputFaceHashes,
                     int                     hashUpperBound)
{
  const char* const   op    = "hollow";
  const TopoDS_Shape& shape = get(kernel, shapeId);
  if (shape.ShapeType() != TopAbs_SOLID || !std::isfinite(thickness) || thickness <= 0.0
      || !std::isfinite(tolerance) || tolerance <= 0.0)
  {
    fail(op, invalid_input, "hollow requires a solid, positive thickness and positive tolerance");
  }
  FaceHistory        history(op, inputFaceHashes, hashUpperBound);
  Lineage            walls(shape);
  const TopoDS_Shape work = withoutVanishing(op, shape, -thickness, walls);
  const TopoDS_Shape inset =
    offsetBody(work, -thickness, tolerance, GeomAbs_Intersection, {}, 0.0, walls);
  if (inset.IsNull())
    fail(op,
         no_solution,
         "the body has no inside left at this wall, or its inside crosses itself; reduce the "
         "thickness");

  Untouched<BRepAlgoAPI_Cut> cut(shape, inset);
  if (!cut.IsDone() || cut.HasErrors() || cut.Shape().IsNull()
      || !BRepCheck_Analyzer(cut.Shape()).IsValid())
  {
    fail(op, no_solution, "the inner cavity could not be subtracted");
  }
  for (TopExp_Explorer faces(shape, TopAbs_FACE); faces.More(); faces.Next())
  {
    if (cut.IsDeleted(faces.Current()) || !cut.Modified(faces.Current()).IsEmpty())
    {
      fail(op, no_solution, "the inner cavity crosses the outer skin", {faces.Current()});
    }
  }
  // A closed cavity adds an inner shell while preserving every original outer face.
  int before = 0, after = 0;
  for (TopExp_Explorer shells(shape, TopAbs_SHELL); shells.More(); shells.Next())
    before++;
  for (TopExp_Explorer shells(cut.Shape(), TopAbs_SHELL); shells.More(); shells.Next())
    after++;
  if (after <= before)
    fail(op, no_solution, "no enclosed cavity was created");
  walls.through(cut);
  walls.within(cut.Shape());
  // Apart from the kernel's own word: the wall is nowhere thinner than asked.
  if (clearance(walls.all(), shape) < thickness - clearanceSlack * tolerance)
    fail(op, no_solution, "the wall comes out thinner than asked; the inside meets itself");
  walls.generated(history);
  return history.finish(store(kernel, cut.Shape()));
}

} // namespace BRepExt
