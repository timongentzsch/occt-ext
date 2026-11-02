// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_OffsetSurface.hxx>
#include <NCollection_HSequence.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <algorithm>
#include <vector>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>

// What is shown in red when a Boolean or a Hollow fails.
// Neither comes out of the failed solve, which reports nothing of where: both are built apart from
// it and say what the operation could not have had, not what it tripped on.

namespace BRepExt
{

// The curves where two bodies meet that do not close: a sheet that stops short of crossing a solid
// meets it in an open section, which cuts nothing. By the Jordan curve theorem on the surface a
// section can only divide a face if it closes on itself or ends on the face's boundary, so these
// are the part of the section that keeps a Boolean from finding a region. What is left of a section
// after joining all its edges; empty where every wire closes.
EvolutionData openSections(OcctKernel& kernel, uint32_t a, uint32_t b, double tolerance)
{
  const char* const op = "openSections";
  if (!std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "tolerance must be positive");
  BRepAlgoAPI_Section section(get(kernel, a), get(kernel, b), false);
  section.ComputePCurveOn1(false);
  section.ComputePCurveOn2(false);
  section.Build();
  if (!section.IsDone())
    fail(op, no_solution, "the two bodies could not be intersected");
  Handle(NCollection_HSequence<TopoDS_Shape>) edges = new NCollection_HSequence<TopoDS_Shape>();
  for (TopExp_Explorer it(section.Shape(), TopAbs_EDGE); it.More(); it.Next())
    edges->Append(it.Current());
  TopoDS_Compound out;
  BRep_Builder    builder;
  builder.MakeCompound(out);
  if (!edges->IsEmpty())
  {
    Handle(NCollection_HSequence<TopoDS_Shape>) wires;
    ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edges, tolerance, false, wires);
    for (int i = 1; i <= wires->Length(); ++i)
    {
      const TopoDS_Wire wire = TopoDS::Wire(wires->Value(i));
      if (!BRep_Tool::IsClosed(wire))
        builder.Add(out, wire);
    }
  }
  return stored(kernel, out);
}

// Each face of a body moved along its normal by `distance`, on its own and untrimmed to the others:
// the shell an offset would have to join. Its surface is the exact parallel surface of the face's,
// over the same parameters; where it passes a focal point it folds, which is what can be seen.
EvolutionData offsetShellFaces(OcctKernel& kernel, uint32_t body, double distance)
{
  const char* const op = "offsetShellFaces";
  if (!std::isfinite(distance))
    fail(op, invalid_input, "the distance has to be a number");
  TopoDS_Compound out;
  BRep_Builder    builder;
  builder.MakeCompound(out);
  for (TopExp_Explorer it(get(kernel, body), TopAbs_FACE); it.More(); it.Next())
  {
    const TopoDS_Face face = TopoDS::Face(it.Current());
    double            u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    // Along the face's outward normal, as the surface's own is reversed with the face.
    const double               along = face.Orientation() == TopAbs_REVERSED ? -distance : distance;
    Handle(Geom_OffsetSurface) surface =
      new Geom_OffsetSurface(BRep_Tool::Surface(face), along, true);
    BRepBuilderAPI_MakeFace made(surface, u0, u1, v0, v1, Precision::Confusion());
    // Where the parallel surface folds there is no face to show.
    if (made.IsDone())
      builder.Add(out, made.Face());
  }
  return stored(kernel, out);
}

// Why a Thicken failed, where that can be told apart from the failed solve: the edges that a move
// of every face by `thickness` shrinks to nothing or turns over. A vertex moves with its faces: to
// the least move that takes it `thickness` along each face's normal there. For flat faces this is
// exact: an edge whose ends pass each other has no place in the offset. On curved faces it uses the
// tangent planes at the vertices, so it is a first-order estimate, and a failure for any other
// reason names no edges. It always throws: the failure, at those edges.
EvolutionData thickenFailure(OcctKernel& kernel, uint32_t sheet, double thickness)
{
  const char* const   op    = "thicken";
  const TopoDS_Shape& shape = get(kernel, sheet);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOf;
  TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_FACE, facesOf);
  // The least-norm move with n_i . p = t, by the normal equations over the normals' span.
  const auto moved = [&](const TopoDS_Vertex& vertex) {
    const gp_Pnt        at = BRep_Tool::Pnt(vertex);
    std::vector<gp_XYZ> normals;
    for (const auto& owner : facesOf.FindFromKey(vertex))
    {
      const TopoDS_Face                   face = TopoDS::Face(owner);
      const Handle(ShapeAnalysis_Surface) on = new ShapeAnalysis_Surface(BRep_Tool::Surface(face));
      const gp_Pnt2d                      uv = on->ValueOfUV(at, Precision::Confusion());
      gp_Pnt                              p;
      gp_Vec                              du, dv;
      BRep_Tool::Surface(face)->D1(uv.X(), uv.Y(), p, du, dv);
      gp_Vec n = du.Crossed(dv);
      if (n.Magnitude() <= gp::Resolution())
        continue;
      n.Normalize();
      if (face.Orientation() == TopAbs_REVERSED)
        n.Reverse();
      if (std::none_of(normals.begin(), normals.end(), [&](const gp_XYZ& o) {
            return o.IsEqual(n.XYZ(), 1e-9);
          }))
        normals.push_back(n.XYZ());
    }
    const int m = int(normals.size());
    if (m == 0)
      return at;
    math_Matrix gram(1, m, 1, m);
    math_Vector right(1, m, thickness), weights(1, m);
    for (int i = 1; i <= m; ++i)
      for (int j = 1; j <= m; ++j)
        gram(i, j) = normals[size_t(i - 1)].Dot(normals[size_t(j - 1)]);
    math_Gauss solve(gram);
    if (!solve.IsDone())
      return at;
    solve.Solve(right, weights);
    gp_XYZ move(0, 0, 0);
    for (int i = 1; i <= m; ++i)
      move += normals[size_t(i - 1)] * weights(i);
    return at.Translated(gp_Vec(move));
  };
  std::vector<TopoDS_Shape>                                     edges;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> all;
  TopExp::MapShapes(shape, TopAbs_EDGE, all);
  for (int i = 1; i <= all.Extent(); ++i)
  {
    const TopoDS_Edge edge = TopoDS::Edge(all(i));
    TopoDS_Vertex     first, last;
    TopExp::Vertices(edge, first, last);
    if (first.IsNull() || last.IsNull() || first.IsSame(last))
      continue;
    const gp_Vec before(BRep_Tool::Pnt(first), BRep_Tool::Pnt(last));
    const gp_Vec after(moved(first), moved(last));
    if (after.Dot(before) <= 0)
      edges.push_back(edge);
  }
  if (edges.empty())
    fail(op, no_solution, "the sheet could not be thickened by this much");
  fail(op,
       no_solution,
       "the sheet could not be thickened by this much: these edges shrink to nothing",
       edges);
}

} // namespace BRepExt
