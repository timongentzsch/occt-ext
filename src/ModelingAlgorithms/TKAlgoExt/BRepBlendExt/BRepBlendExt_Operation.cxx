// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The blend as an operation on a stored body: its edges are named by handle, the request is made,
// OpenCascade's own failures are reported by name, and the result is accepted only if it is a
// valid, closed body whose curvature-continuous seams hold their claim.

#include <BRepToolsExt_Session.hxx>
#include <BRepFilletExt_Blend.hxx>
#include <BRepBlendExt_Internal.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRep_Tool.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

namespace BRepExt
{

namespace blending
{
static int freeEdges(const TopoDS_Shape& shape)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    faces;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, faces);
  int count = 0;
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    const TopoDS_Edge edge = TopoDS::Edge(faces.FindKey(i));
    const auto&       on   = faces.FindFromIndex(i);
    count += on.Extent() == 1 && !BRep_Tool::Degenerated(edge)
             && !BRep_Tool::IsClosed(edge, TopoDS::Face(on.First()));
  }
  return count;
}

Blended blendStored(const OcctKernel&            kernel,
                    const char*                  op,
                    uint32_t                     shape,
                    const std::vector<uint32_t>& edges,
                    double                       first,
                    double                       second,
                    int                          law,
                    double                       rho,
                    uint32_t                     reference,
                    const std::vector<int>&      counts,
                    const std::vector<double>&   sizes,
                    double                       angle,
                    bool                         deferredFit)
{
  if (law < 0 || law > 4)
    fail(op, invalid_input, "unknown section law");
  if (!counts.empty() && counts.size() != edges.size())
    fail(op, invalid_input, "one count of size places for each edge");
  Request request;
  request.operation = op;
  request.fitWork   = deferredFit ? FitWork::Deferred : FitWork::Bounded;
  request.body      = get(kernel, shape);
  for (const auto& edge : subShapes(kernel, op, request.body, edges, TopAbs_EDGE))
    request.edges.push_back(edge);
  request.first  = first;
  request.second = second;
  request.law    = static_cast<Law>(law);
  request.rho    = rho;
  request.angle  = angle;
  if (reference)
    request.reference = subShapes(kernel, op, request.body, {reference}, TopAbs_FACE)[0];
  size_t next = 0;
  for (const int count : counts)
  {
    auto& points = request.radii.emplace_back();
    for (int k = 0; k < count; ++k, next += 2)
    {
      if (next + 1 >= sizes.size())
        fail(op, invalid_input, "fewer size places than counted");
      points.emplace_back(sizes[next], sizes[next + 1]);
    }
  }
  if (next != sizes.size())
    fail(op, invalid_input, "more size places than counted");
  Blended blended;
  try
  {
    blended = blend(request);
  }
  catch (const Standard_Failure& error)
  {
    // OpenCascade's own failures carry its type, often no text: name both.
    fail(op,
         no_solution,
         std::string("OpenCascade could not build the blend (") + error.ExceptionType()
           + (*error.what() ? std::string(": ") + error.what() : std::string()) + ")");
  }
  // The blend's size: its largest, where it changes along the edges.
  double size = std::max(first, second);
  for (const auto& points : request.radii)
    for (const auto& [place, at] : points)
      size = std::max(size, std::max(first, second) * at / first);
  accepted(op, request.body, blended.shape, size, &blended);
  return blended;
}

void accepted(const char*         op,
              const TopoDS_Shape& input,
              const TopoDS_Shape& result,
              double              size,
              const Blended*      made)
{
  if (const auto what = invalidity(result, input); !what.empty())
    fail(op,
         no_solution,
         "the blended body is not valid (shape kind and OpenCascade status:" + what + ")");
  // OpenCascade's check passes a shell with an edge on one face: none that the input did not have.
  if (freeEdges(result) > freeEdges(input))
    fail(op, no_solution, "the blended body is not closed: an edge of it lies on one face");
  // Every edge that claims curvature continuity holds it (Attach sets the claims), read off the
  // delivered surfaces with the measure of Measure Continuity: three places inside every span of
  // the edge, so that a short span is seen, and the G2 oracle's seventeen on a curve with none.
  // Not an edge kept as it was, whose claim an earlier blend of another size made.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> kept;
  if (made)
    for (const auto& [before, after] : made->unchanged)
      kept.Add(after);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    faces;
  TopExp::MapShapesAndAncestors(result, TopAbs_EDGE, TopAbs_FACE, faces);
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    const TopoDS_Edge edge = TopoDS::Edge(faces.FindKey(i));
    const auto&       on   = faces.FindFromIndex(i);
    if (on.Extent() != 2 || BRep_Tool::Degenerated(edge) || kept.Contains(edge))
      continue;
    const TopoDS_Face a = TopoDS::Face(on.First()), b = TopoDS::Face(on.Last());
    if (a.IsSame(b) || !BRep_Tool::HasContinuity(edge, a, b)
        || BRep_Tool::Continuity(edge, a, b) < GeomAbs_G2)
      continue;
    const BRepAdaptor_Curve    curve(edge);
    const BoundaryJet          ja(op, edge, a), jb(op, edge, b);
    const int                  spans = curve.NbIntervals(GeomAbs_CN);
    NCollection_Array1<double> ends(1, spans + 1);
    curve.Intervals(ends, GeomAbs_CN);
    const int samples = spans == 1 ? 17 : 3;
    for (int k = 0; k < spans * samples; ++k)
    {
      const int    span = 1 + k / samples;
      const double u = ends(span) + (ends(span + 1) - ends(span)) * (k % samples + 0.5) / samples;
      gp_Pnt       p;
      gp_Vec       t;
      curve.D1(u, p, t);
      if (!(t.Magnitude() > detail::limit::degenerate))
        continue;
      const double bent = bendDifference(ja.at(u), jb.at(u), t.Normalized()).first;
      if (bent * size > detail::limit::seamBend
          && !detail::unresolved(bent, t.Magnitude() * (ends(span + 1) - ends(span)), size))
        fail(op,
             no_solution,
             "a seam of the G2 blend is not curvature continuous",
             {edge},
             FailureKind::Curvature);
    }
  }
}
} // namespace blending

} // namespace BRepExt
