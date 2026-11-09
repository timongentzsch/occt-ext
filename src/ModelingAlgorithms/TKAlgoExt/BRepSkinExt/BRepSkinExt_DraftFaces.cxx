// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
// draftFaces: taper faces for mould release.
//
// BRepOffsetAPI_DraftAngle is OCCT's only taper algorithm. It inclines planar, cylindrical and
// conical faces about their intersection with a neutral plane and re-intersects the neighbours.
// One maker takes any number of faces, each validated on its own by AddDone.

#include <BRepSkinExt_DraftIsocline.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepOffsetExt_Solid.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <gp_Cylinder.hxx>
#include <BRepOffsetAPI_DraftAngle.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <optional>
#include <cmath>

namespace BRepExt
{

namespace
{

const char* const op = "draftFaces";

//! The tolerance that the draft of curved walls works to: the operation is given none, and
//! this is the one that documents are made with.
constexpr double draftTolerance = BRepExt::tolerance::sewing;

gp_Dir direction(Vec3 v, const char* name)
{
  if (!finite(v))
  {
    fail(op, invalid_input, std::string("the ") + name + " must be finite");
  }
  const gp_Vec vec(v.x, v.y, v.z);
  if (vec.Magnitude() <= gp::Resolution())
  {
    fail(op, invalid_input, std::string("the ") + name + " is a zero vector");
  }
  return gp_Dir(vec);
}

//! The faces to draft split into those OpenCascade's draft takes (`taken`: planes, cylinders and
//! cones) and the walls of any other kind, drafted first by the ruled isocline surface. A flat or
//! round wall that goes over smoothly into such a wall is drafted with it, the same way: the two
//! share the line they meet in, before and after.
void splitWalls(const TopoDS_Shape&              shape,
                const std::vector<TopoDS_Shape>& faces,
                std::vector<TopoDS_Shape>&       walls,
                std::vector<TopoDS_Shape>&       taken)
{
  for (const TopoDS_Shape& face : faces)
  {
    const GeomAbs_SurfaceType type = BRepAdaptor_Surface(TopoDS::Face(face)).GetType();
    (type == GeomAbs_Plane || type == GeomAbs_Cylinder || type == GeomAbs_Cone ? taken : walls)
      .push_back(face);
  }
  // A flat or round wall that goes over smoothly into such a wall is drafted with it, the
  // same way: the two share the line they meet in, before and after.
  if (!walls.empty())
  {
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>
      facesOf;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, facesOf);
    for (size_t next = 0; next < walls.size(); next++)
      for (TopExp_Explorer edges(walls[next], TopAbs_EDGE); edges.More(); edges.Next())
        for (const TopoDS_Shape& other : facesOf.FindFromKey(edges.Current()))
        {
          const auto found =
            std::find_if(taken.begin(), taken.end(), [&](const TopoDS_Shape& face) {
              return face.IsSame(other);
            });
          if (found == taken.end()
              || BRepAdaptor_Surface(TopoDS::Face(other)).GetType() == GeomAbs_Cone
              || !tangentAlong(TopoDS::Edge(edges.Current()),
                               TopoDS::Face(walls[next]),
                               TopoDS::Face(other)))
            continue;
          walls.push_back(*found);
          taken.erase(found);
        }
  }
}

//! The walls drafted by the isocline surface, and the faces OpenCascade takes drafted beside
//! them.
EvolutionData draftCurvedWalls(OcctKernel&                      kernel,
                               const TopoDS_Shape&              shape,
                               const std::vector<TopoDS_Shape>& walls,
                               const std::vector<TopoDS_Shape>& taken,
                               double                           angle,
                               const gp_Dir&                    pull,
                               const gp_Pln&                    neutral,
                               FaceHistory&                     history)
{
  const TopoDS_Shape solid = singleSolid(shape);
  if (solid.IsNull())
    fail(op, unsupported_geometry, "a curved face is drafted on a solid only", walls);
  Lineage      lineage(solid);
  TopoDS_Shape result = draftWalls(solid, walls, angle, pull, neutral, draftTolerance, lineage);
  if (!taken.empty())
  {
    BRepOffsetAPI_DraftAngle rest(result);
    for (const TopoDS_Shape& face : taken)
    {
      const auto& now = lineage.of(face);
      if (now.size() != 1)
        fail(op,
             no_solution,
             "a face to draft did not come through the draft of the walls beside it as one face",
             {face});
      rest.Add(TopoDS::Face(now[0]), pull, angle, neutral);
      if (!rest.AddDone())
        fail(op,
             unsupported_geometry,
             "a face cannot be drafted with this pull direction and neutral plane",
             {face});
    }
    rest.Build();
    if (!rest.IsDone() || !BRepCheck_Analyzer(rest.Shape()).IsValid())
      fail(op,
           no_solution,
           "the flat, cylindrical and conical faces could not be drafted beside the drafted walls",
           taken);
    lineage.throughBoth(rest);
    result = rest.Shape();
    lineage.within(result);
  }
  lineage.modified(history);
  return history.finish(store(kernel, result));
}

//! OCCT can fail to propagate a planar wall's draft through tangent round
//! walls. The ruled isocline construction can draft the complete tangent chain;
//! fails with why `face` (`which`) could not be drafted otherwise.
EvolutionData tangentChainDraft(OcctKernel&                      kernel,
                                const TopoDS_Shape&              shape,
                                const std::vector<TopoDS_Shape>& faces,
                                const TopoDS_Face&               face,
                                const std::string&               which,
                                double                           angle,
                                const gp_Dir&                    pull,
                                const gp_Pln&                    neutral,
                                FaceHistory&                     history)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    adjacent;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacent);
  std::vector<TopoDS_Shape> chain = faces;
  for (size_t next = 0; next < chain.size(); next++)
    for (TopExp_Explorer edges(chain[next], TopAbs_EDGE); edges.More(); edges.Next())
      for (const TopoDS_Shape& other : adjacent.FindFromKey(edges.Current()))
      {
        if (std::any_of(chain.begin(), chain.end(), [&](const TopoDS_Shape& have) {
              return have.IsSame(other);
            }))
          continue;
        const BRepAdaptor_Surface support(TopoDS::Face(other));
        const bool                wall =
          (support.GetType() == GeomAbs_Plane
           && std::abs(support.Plane().Axis().Direction().Dot(pull)) < alongSlack)
          || (support.GetType() == GeomAbs_Cylinder
              && std::abs(support.Cylinder().Axis().Direction().Dot(pull)) > 1.0 - alongSlack);
        if (wall
            && tangentAlong(TopoDS::Edge(edges.Current()),
                            TopoDS::Face(chain[next]),
                            TopoDS::Face(other)))
          chain.push_back(other);
      }
  const TopoDS_Shape        solid = singleSolid(shape);
  const BRepAdaptor_Surface selected(face);
  if (selected.GetType() == GeomAbs_Plane
      && std::abs(selected.Plane().Axis().Direction().Dot(pull)) < alongSlack)
  {
    for (const TopoDS_Shape& wall : chain)
      for (TopExp_Explorer edges(wall, TopAbs_EDGE); edges.More(); edges.Next())
        for (const TopoDS_Shape& other : adjacent.FindFromKey(edges.Current()))
        {
          if (std::any_of(chain.begin(), chain.end(), [&](const TopoDS_Shape& have) {
                return have.IsSame(other);
              }))
            continue;
          if (tangentAlong(TopoDS::Edge(edges.Current()), TopoDS::Face(wall), TopoDS::Face(other)))
            fail(op,
                 unsupported_geometry,
                 "a tangent rounded cap does not run along the pull direction; this draft "
                 "cannot rebuild its junction with the wall",
                 {wall, other});
        }
  }
  if (chain.size() > faces.size() && !solid.IsNull())
  {
    Lineage            lineage(solid);
    const TopoDS_Shape result =
      draftWalls(solid, chain, angle, pull, neutral, draftTolerance, lineage);
    lineage.modified(history);
    return history.finish(store(kernel, result));
  }
  const bool parallel = BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane
                        && std::abs(BRepAdaptor_Surface(face).Plane().Axis().Direction().Dot(
                             neutral.Axis().Direction()))
                             > 1.0 - alongSlack;
  fail(op,
       unsupported_geometry,
       which
         + (parallel ? " is parallel to the neutral plane and cannot be inclined about an "
                       "intersection with it"
                     : " could not be drafted: OpenCascade could not recompute the wall and "
                       "its tangent neighbours"),
       {face});
}

//! DraftAngle preserves adjacency and can leave an invalid wall at breakout.
//! Recut an expanding, complete cylindrical cavity by its exact conical support.
//! Boundaries must be perpendicular planar caps; shrinking holes and oblique
//! pull directions need a different reconstruction and never enter this path.
std::optional<EvolutionData> conicalRecut(OcctKernel&         kernel,
                                          const TopoDS_Shape& shape,
                                          const TopoDS_Face&  wall,
                                          double              angle,
                                          const gp_Dir&       pull,
                                          const gp_Pln&       neutral,
                                          FaceHistory&        history)
{
  BRepAdaptor_Surface s(wall);
  if (s.GetType() == GeomAbs_Cylinder && wall.Orientation() == TopAbs_REVERSED
      && std::abs(s.LastUParameter() - s.FirstUParameter() - 2 * M_PI)
           < BRepExt::tolerance::angularResolution)
  {
    const auto   cylinder  = s.Cylinder();
    const auto   axis      = cylinder.Axis().Direction();
    const double alignment = axis.Dot(pull);
    bool         eligible =
      std::abs(alignment) > 1 - BRepExt::tolerance::angle
      && std::abs(axis.Dot(neutral.Axis().Direction())) > 1 - BRepExt::tolerance::angle;
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>
      adj;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adj);
    for (TopExp_Explorer e(wall, TopAbs_EDGE); e.More(); e.Next())
      for (const auto& f : adj.FindFromKey(e.Current()))
        if (!f.IsSame(wall))
        {
          BRepAdaptor_Surface other(TopoDS::Face(f));
          if (other.GetType() != GeomAbs_Plane
              || std::abs(other.Plane().Axis().Direction().Dot(axis))
                   < 1 - BRepExt::tolerance::angle)
            eligible = false;
        }
    const auto   start  = cylinder.Location().Translated(gp_Vec(axis) * s.FirstVParameter());
    const double height = s.LastVParameter() - s.FirstVParameter();
    const double axial  = gp_Vec(neutral.Location(), start).Dot(gp_Vec(pull));
    const double r0     = cylinder.Radius() + axial * std::tan(angle);
    const double r1     = r0 + height * alignment * std::tan(angle);
    if (eligible && r0 >= cylinder.Radius() - Precision::Confusion()
        && r1 >= cylinder.Radius() - Precision::Confusion() && height > Precision::Confusion())
    {
      const auto tool = BRepPrimAPI_MakeCone(gp_Ax2(start, axis), r0, r1, height).Shape();
      Untouched<BRepAlgoAPI_Cut> cut(shape, tool);
      if (cut.IsDone() && BRepCheck_Analyzer(cut.Shape()).IsValid())
      {
        history.collect(cut, shape);
        // The new conical faces descend from this wall, through the tool.
        NCollection_List<TopoDS_Shape> successors;
        for (TopExp_Explorer f(tool, TopAbs_FACE); f.More(); f.Next())
        {
          if (BRepAdaptor_Surface(TopoDS::Face(f.Current())).GetType() != GeomAbs_Cone)
            continue;
          for (const auto& next : cut.Modified(f.Current()))
            if (next.ShapeType() == TopAbs_FACE)
              successors.Append(next);
        }
        if (!successors.IsEmpty())
          history.modified(wall, successors);
        return history.finish(store(kernel, cut.Shape()));
      }
    }
  }
  return std::nullopt;
}

} // namespace

EvolutionData draftFaces(OcctKernel&                  kernel,
                         uint32_t                     shapeId,
                         const std::vector<uint32_t>& faceIds,
                         double                       angle,
                         Vec3                         pullDirection,
                         Vec3                         planeOrigin,
                         Vec3                         planeNormal,
                         const std::vector<int>&      inputFaceHashes,
                         int                          hashUpperBound)
{
  const TopoDS_Shape& shape = get(kernel, shapeId);
  if (!finite(planeOrigin))
  {
    fail(op, invalid_input, "the neutral plane origin must be finite");
  }
  const gp_Dir pull = direction(pullDirection, "pull direction");
  const gp_Pln neutral(gp_Pnt(planeOrigin.x, planeOrigin.y, planeOrigin.z),
                       direction(planeNormal, "neutral plane normal"));
  if (!std::isfinite(angle) || std::abs(angle) <= Precision::Angular()
      || std::abs(angle) >= M_PI / 2)
  {
    fail(op, invalid_input, "the angle must be non-zero and below 90 degrees");
  }
  FaceHistory history(op, inputFaceHashes, hashUpperBound);

  const std::vector<TopoDS_Shape> faces = subShapes(kernel, op, shape, faceIds, TopAbs_FACE);
  // OpenCascade's draft takes planes, cylinders and cones, and ignores other faces without
  // reporting it. Walls of any other kind are drafted first, by the ruled isocline surface.
  std::vector<TopoDS_Shape> walls, taken;
  splitWalls(shape, faces, walls, taken);
  if (!walls.empty())
    return draftCurvedWalls(kernel, shape, walls, taken, angle, pull, neutral, history);
  BRepOffsetAPI_DraftAngle maker(shape);
  for (size_t i = 0; i < faces.size(); i++)
  {
    const TopoDS_Face& face = TopoDS::Face(faces[i]);
    const std::string  which =
      "face " + std::to_string(i + 1) + " of " + std::to_string(faces.size());
    maker.Add(face, pull, angle, neutral);
    if (!maker.AddDone())
      return tangentChainDraft(kernel, shape, faces, face, which, angle, pull, neutral, history);
  }
  maker.Build();
  if ((!maker.IsDone() || !BRepCheck_Analyzer(maker.Shape()).IsValid()) && faces.size() == 1)
    if (const auto recut =
          conicalRecut(kernel, shape, TopoDS::Face(faces[0]), angle, pull, neutral, history))
      return *recut;
  if (!maker.IsDone())
  {
    fail(op, no_solution, "the drafted faces could not be re-intersected with their neighbours");
  }
  // Draft_Modification reports "done" for a taper that runs out of its neighbours.
  const TopoDS_Shape& result = maker.Shape();
  if (!BRepCheck_Analyzer(result).IsValid())
  {
    fail(op, no_solution, "the taper runs out of the neighbouring faces; reduce the angle");
  }
  // The maker files a drafted face's successor under Generated and a re-trimmed neighbour's
  // under Modified. Both are the same face afterwards, so both are reported as modified.
  for (TopExp_Explorer all(shape, TopAbs_FACE); all.More(); all.Next())
  {
    const TopoDS_Shape& face = all.Current();
    if (!history.tracks(face))
      continue;
    NCollection_List<TopoDS_Shape> successors = maker.Modified(face);
    for (const TopoDS_Shape& successor : maker.Generated(face))
      successors.Append(successor);
    history.modified(face, successors);
  }
  return history.finish(store(kernel, result));
}

} // namespace BRepExt
