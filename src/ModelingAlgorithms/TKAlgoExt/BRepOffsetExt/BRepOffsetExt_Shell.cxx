// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// A solid hollowed to a wall, with faces of it removed to open it.
//
// Where no face vanishes OpenCascade's own shell (BRepOffsetAPI_MakeThickSolid) is tried
// first. It cannot remove a blend whose radius the wall reaches, nor open a face that is
// tangent to its neighbours. Then, and where it returns something that does not hold:
//
//   1. the faces that the wall makes vanish are removed from a working copy
//   (BRepOffsetExt_Solid.cxx);
//   2. the copy is offset inwards, with the faces to open moved a little outwards instead, so
//      that the walls beside an opening run on through it; where an opening is tangent to
//      its neighbours they cannot run on, and the copy is offset as a whole and a prism of
//      the opening's own image stands in for it;
//   3. that body is subtracted from the ORIGINAL, whose outer blends therefore stay, and since
//      it stands out of the faces to open, what is left of those is the rim of the opening.
//
// Local (hollow every face but those that border
// no opening): only the faces round the openings get a wall. The others are removed from the
// working copy first and the faces round the openings extended to close it, so that the inside
// is that of the walls round the openings alone, subtracted from the original: what lies behind
// the other faces stays solid, a pocket there that reaches the inside opens into it. Where the
// faces round the openings cannot close it, nothing is made.

#include <BRepOffsetExt_Solid.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepFilletExt_Untouched.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <cmath>
#include <optional>

namespace BRepExt
{

namespace
{

const char* const op = "shell";

//! OpenCascade's shell, if what it returns holds.
std::optional<EvolutionData> quickShell(OcctKernel&                      kernel,
                                        const TopoDS_Shape&              source,
                                        const std::vector<TopoDS_Shape>& removed,
                                        double                           thickness,
                                        double                           tolerance,
                                        int                              join,
                                        const Lineage&                   divided,
                                        FaceHistory&                     history)
{
  NCollection_List<TopoDS_Shape> open;
  for (const TopoDS_Shape& face : removed)
    open.Append(face);
  try
  {
    BRepOffsetAPI_MakeThickSolid maker;
    maker.MakeThickSolidByJoin(source,
                               open,
                               -thickness,
                               tolerance,
                               BRepOffset_Skin,
                               false,
                               false,
                               static_cast<GeomAbs_JoinType>(join),
                               false,
                               progress(true)->Start());
    checkCancelled(op);
    if (!maker.IsDone() || maker.Shape().IsNull() || !BRepCheck_Analyzer(maker.Shape()).IsValid())
      return {};
    // A wall below what the offset resolves comes back as the input.
    const double before = estimatedVolume(source), after = estimatedVolume(maker.Shape());
    if (std::abs(after) <= 1e-15 || std::abs(after - before) <= std::abs(before) * 1e-8)
      return {};
    // An inward wall takes material away and an outward one is all new.
    if (thickness > 0 ? after >= before : after <= 0)
      return {};
    Lineage walls(source);
    walls.images(maker);
    walls.within(maker.Shape());
    if (!clearance(walls.all(removed), source, removed, std::abs(thickness),
                   std::abs(thickness) - clearanceSlack * tolerance)
           .satisfies(std::abs(thickness) - clearanceSlack * tolerance))
      return {};
    Lineage changed(source);
    changed.through(maker);
    changed.within(maker.Shape());
    Lineage outer = divided, inner = divided;
    outer.then(changed);
    inner.then(walls);
    outer.modified(history);
    inner.generated(history);
    return history.finish(store(kernel, maker.Shape()));
  }
  catch (const Standard_Failure&)
  {
    checkCancelled(op);
    return {};
  }
}

bool tangentToANeighbour(const TopoDS_Shape& body, const TopoDS_Shape& face)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOf;
  TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOf);
  for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
    for (const TopoDS_Shape& other : facesOf.FindFromKey(edges.Current()))
      if (!other.IsSame(face)
          && tangentAlong(TopoDS::Edge(edges.Current()), TopoDS::Face(face), TopoDS::Face(other)))
        return true;
  return false;
}

//! The faces of `body` that border none of `openings` and are none of them.
std::vector<TopoDS_Shape> awayFrom(const TopoDS_Shape&              body,
                                   const std::vector<TopoDS_Shape>& openings)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOf;
  TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOf);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> near;
  for (const TopoDS_Shape& opening : openings)
    for (TopExp_Explorer edges(opening, TopAbs_EDGE); edges.More(); edges.Next())
      for (const TopoDS_Shape& face : facesOf.FindFromKey(edges.Current()))
        near.Add(face);
  std::vector<TopoDS_Shape>                                     away;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> faces;
  TopExp::MapShapes(body, TopAbs_FACE, faces);
  for (int i = 1; i <= faces.Extent(); i++)
    if (!near.Contains(faces(i)))
      away.push_back(faces(i));
  return away;
}

} // namespace

EvolutionData shell(OcctKernel&                  kernel,
                    uint32_t                     id,
                    const std::vector<uint32_t>& ids,
                    double                       thickness,
                    double                       tolerance,
                    int                          join,
                    const std::vector<int>&      hashes,
                    int                          bound,
                    bool                         local)
{
  const TopoDS_Shape given = singleSolid(get(kernel, id));
  if (given.IsNull() || !std::isfinite(thickness) || thickness == 0 || !std::isfinite(tolerance)
      || tolerance <= 0 || join < 0 || join > 2)
    fail(op,
         invalid_input,
         "one solid, a thickness that is not zero, a positive tolerance and a supported join are "
         "required");
  // The measured wall can be clearanceSlack * tolerance thinner than requested. When that
  // allowance consumes the whole wall, the check cannot establish that a wall remains.
  if (std::abs(thickness) <= clearanceSlack * tolerance)
    fail(op,
         no_solution,
         "the wall is no larger than the clearance allowance; increase the wall thickness or lower "
         "the offset tolerance");
  FaceHistory history(op, hashes, bound);
  // Divided along any C0 surface's knots first (dividedAtC0); the faces to open go with it.
  Lineage                   divided(given);
  const TopoDS_Shape        source = dividedAtC0(given, divided);
  std::vector<TopoDS_Shape> removed;
  for (const TopoDS_Shape& face : subShapes(kernel, op, given, ids, TopAbs_FACE))
    for (const TopoDS_Shape& part : divided.of(face))
      removed.push_back(part);

  Lineage inner(source);
  // With Local, the body whose inside is hollowed: the source without the faces left alone.
  TopoDS_Shape              closed = source;
  std::vector<TopoDS_Shape> alone, openedInClosed = removed;
  if (local)
  {
    if (removed.empty())
      fail(op,
           invalid_input,
           "Local hollows round the faces to open; pick the faces to open, or turn Local off");
    alone = awayFrom(source, removed);
    if (alone.empty())
      fail(op,
           no_solution,
           "every face borders a face to open, so Local leaves none alone; turn Local off to "
           "hollow the whole body");
    BRepAlgoAPI_Defeaturing defeaturing;
    defeaturing.SetShape(source);
    for (const TopoDS_Shape& face : alone)
      defeaturing.AddFaceToRemove(face);
    defeaturing.SetToFillHistory(true);
    defeaturing.Build(progress(true)->Start());
    checkCancelled(op);
    closed    = defeaturing.IsDone() && !defeaturing.HasErrors() ? singleSolid(defeaturing.Shape())
                                                                 : TopoDS_Shape();
    bool gone = !closed.IsNull();
    for (const TopoDS_Shape& face : alone)
      gone = gone && defeaturing.IsDeleted(face);
    if (!gone || !BRepCheck_Analyzer(closed).IsValid())
      fail(op,
           no_solution,
           "the faces round the openings do not close a hollow of their own: it would run into the "
           "faces "
           "that Local leaves alone",
           alone);
    inner.through(defeaturing);
    openedInClosed.clear();
    for (const TopoDS_Shape& face : removed)
      for (const TopoDS_Shape& image : inner.of(face))
        openedInClosed.push_back(image);
  }
  // Fails where a face is too tight for the wall in a part of it only.
  const TopoDS_Shape work            = withoutVanishing(op, closed, -thickness, inner);
  const bool         nothingVanishes = work.IsSame(closed);
  // It cannot open a face that is tangent to its neighbours (it hands back the input), so it
  // is not asked to.
  const bool openingTangent =
    std::any_of(removed.begin(), removed.end(), [&](const TopoDS_Shape& face) {
      return tangentToANeighbour(source, face);
    });
  if (nothingVanishes && !local && !openingTangent)
    if (auto quick =
          quickShell(kernel, source, removed, thickness, tolerance, join, divided, history))
      return *quick;
  if (thickness < 0)
    fail(op,
         no_solution,
         nothingVanishes ? "the wall outside the body could not be built"
                         : "a wall outside the body that makes faces vanish is not built",
         removed);

  // The faces to open, as the working copy has them.
  std::vector<TopoDS_Shape> openings;
  bool                      tangent = false;
  for (const TopoDS_Shape& face : removed)
  {
    const auto& now = inner.of(face);
    if (now.size() != 1)
      fail(op, unsupported_geometry, "a face to open is itself gone at this wall", {face});
    openings.push_back(now[0]);
    tangent = tangent || tangentToANeighbour(work, now[0]);
  }
  const double overshoot = std::max(openingOvershoot * thickness, 1e3 * tolerance);

  TopoDS_Shape tool;
  Lineage      walls = inner;
  // The whole inside offset, each opening's face then carried out through the wall as a prism.
  const auto wholeInside = [&] {
    walls = inner;
    tool  = offsetBody(work, -thickness, tolerance, join, {}, 0.0, walls);
    if (tool.IsNull())
      fail(op,
           no_solution,
           "the body has no inside left at this wall, or its inside crosses itself");
    for (size_t i = 0; i < removed.size(); i++)
    {
      const auto& images = walls.of(removed[i]);
      if (images.size() != 1
          || BRepAdaptor_Surface(TopoDS::Face(images[0])).GetType() != GeomAbs_Plane)
        fail(op,
             unsupported_geometry,
             "a curved face that is tangent to its neighbours cannot be opened: the walls "
             "beside it do not run on through it",
             {removed[i]});
      const TopoDS_Face image  = TopoDS::Face(images[0]);
      gp_Dir            normal = BRepAdaptor_Surface(image).Plane().Axis().Direction();
      if (BRepAdaptor_Surface(image).Plane().Direct() == (image.Orientation() == TopAbs_REVERSED))
        normal.Reverse();
      const TopoDS_Shape prism =
        BRepPrimAPI_MakePrism(image, gp_Vec(normal) * (thickness + overshoot)).Shape();
      Untouched<BRepAlgoAPI_Fuse> joined(tool, prism);
      const TopoDS_Shape          both =
        joined.IsDone() && !joined.HasErrors() ? singleSolid(joined.Shape()) : TopoDS_Shape();
      if (both.IsNull() || !BRepCheck_Analyzer(both).IsValid())
        fail(op, no_solution, "the opening could not be joined to the inside", {removed[i]});
      walls.through(joined);
      tool = both;
    }
  };
  // Pushing the openings out by joins of intersection keeps OpenCascade for an hour, past any
  // cancelling, over a surface only C1: its offset is only C0, and every edge's tolerance reads its
  // higher derivatives (BRepOffset_MakeOffset's UpdateTolerance). Over one, the whole inside first,
  // and pushing only for a body as it came (a body divided at C0 has kept it longest).
  bool onlyC1 = false;
  for (TopExp_Explorer f(work, TopAbs_FACE); f.More() && !onlyC1; f.Next())
    onlyC1 = BRep_Tool::Surface(TopoDS::Face(f.Current()))->Continuity() < GeomAbs_C2;
  const auto pushed = [&] {
    walls = inner;
    tool  = tangent ? TopoDS_Shape()
                    : offsetBody(work, -thickness, tolerance, join, openings, overshoot, walls);
  };
  if (onlyC1)
  {
    try
    {
      wholeInside();
    }
    catch (const Failure&)
    {
      checkCancelled(op);
      if (!source.IsSame(given))
        throw;
      pushed();
      if (tool.IsNull())
        throw;
    }
  }
  else
  {
    pushed();
    if (tool.IsNull())
      wholeInside();
  }

  Untouched<BRepAlgoAPI_Cut> cut(source, tool);
  const TopoDS_Shape         result =
    cut.IsDone() && !cut.HasErrors() ? singleSolid(cut.Shape()) : TopoDS_Shape();
  if (result.IsNull() || !BRepCheck_Analyzer(result).IsValid())
    fail(op, no_solution, "the inside could not be subtracted, or the wall falls apart", removed);
  Lineage outer(source);
  outer.through(cut);
  outer.within(result);
  walls.through(cut);
  walls.within(result);
  // An opening that the inside did not reach leaves its face whole.
  for (const TopoDS_Shape& face : removed)
    if (outer.of(face).size() == 1 && outer.of(face)[0].IsSame(face))
      fail(op, no_solution, "the inside does not reach a face that is to be opened", {face});
  // Apart from the kernel's own word: the wall is nowhere thinner than asked. With Local, the
  // wall is that of the faces round the openings: what lies behind the faces left alone is
  // the body's own, which may come nearer or meet the inside.
  // Opening descendants form the rim, not an inward offset wall. They intentionally
  // meet the retained source skin at zero distance. Keep them in the public history,
  // but measure the inward descendants of every retained source face.
  const auto measured=clearance(walls.all(removed), closed, openedInClosed, thickness,
                                thickness-clearanceSlack*tolerance);
  if (!measured.satisfies(thickness-clearanceSlack*tolerance))
    fail(op, no_solution, "the wall clearance is unresolved or below the requested thickness"
         "; minimum="+std::to_string(measured.minimum)+"; covered="+std::to_string(measured.covered)
         +"/"+std::to_string(measured.expected)+"; "+measured.issue);
  // Compose through the initial division: public history is keyed by the caller's faces.
  Lineage originalOuter = divided, originalWalls = divided;
  originalOuter.then(outer);
  originalWalls.then(walls);
  originalOuter.modified(history);
  originalWalls.generated(history);
  return history.finish(store(kernel, result));
}

} // namespace BRepExt
