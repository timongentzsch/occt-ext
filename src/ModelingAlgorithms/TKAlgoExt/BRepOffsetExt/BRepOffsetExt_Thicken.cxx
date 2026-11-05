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

// Sheet thickening shares offset's C0 division and cancellation, including its
// edge-tolerance loop. Dividing changes topology, not the support geometry.
#include <BRepToolsExt_Operations.hxx>
#include <BRepOffsetExt_Parts.hxx>
#include <BRepOffsetExt_Focal.hxx>
#include <TopoDS.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepOffset_MakeOffset.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <cmath>

namespace BRepExt
{
EvolutionData thicken(OcctKernel&             kernel,
                      uint32_t                id,
                      double                  thickness,
                      double                  tolerance,
                      const std::vector<int>& hashes,
                      int                     bound)
{
  const char*        op    = "thicken";
  const TopoDS_Shape given = get(kernel, id);
  if (!std::isfinite(thickness) || thickness == 0 || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "a finite nonzero thickness and positive tolerance are required");
  if (given.ShapeType() == TopAbs_FACE || given.ShapeType() == TopAbs_SHELL)
    for (TopExp_Explorer faces(given, TopAbs_FACE); faces.More(); faces.Next())
      if (focalSweep(TopoDS::Face(faces.Current()), thickness) == FocalSweep::Crosses)
        fail(op,
             no_solution,
             "the thickness passes a principal radius of curvature: the normal sweep folds",
             {faces.Current()},
             FailureKind::Curvature);
  FaceHistory         history(op, hashes, bound);
  BRepBuilderAPI_Copy copy(given, true, false);
  Lineage             divided(given);
  divided.through(copy);
  const TopoDS_Shape work = dividedAtC0(copy.Shape(), divided);
  if (work.ShapeType() == TopAbs_FACE || work.ShapeType() == TopAbs_SHELL)
  {
    BRepOffset_MakeOffset maker;
    maker.Initialize(work, thickness, tolerance, BRepOffset_Skin, false, false, GeomAbs_Arc, true);
    maker.MakeOffsetShape(progress(true)->Start());
    checkCancelled(op);
    if (!maker.IsDone() || maker.Shape().IsNull())
      fail(op, no_solution, "the divided sheet could not be thickened");
    Lineage changed(work), walls(work);
    changed.through(maker);
    walls.images(maker);
    changed.within(maker.Shape());
    walls.within(maker.Shape());
    Lineage outer = divided, inner = divided;
    outer.then(changed);
    inner.then(walls);
    outer.modified(history);
    inner.generated(history);
    return history.finish(store(kernel, maker.Shape()));
  }
  NCollection_List<TopoDS_Shape> empty;
  BRepOffsetAPI_MakeThickSolid   maker;
  maker.MakeThickSolidByJoin(work,
                             empty,
                             thickness,
                             tolerance,
                             BRepOffset_Skin,
                             false,
                             false,
                             GeomAbs_Arc,
                             false,
                             progress(true)->Start());
  checkCancelled(op);
  if (!maker.IsDone() || maker.Shape().IsNull())
    fail(op, no_solution, "the body could not be thickened");
  Lineage changed(work);
  changed.throughBoth(maker);
  changed.within(maker.Shape());
  divided.then(changed);
  divided.modified(history);
  return history.finish(store(kernel, maker.Shape()));
}
} // namespace BRepExt
