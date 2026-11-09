// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepBooleanExt_Accept.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_Timing.hxx>

#include <BOPAlgo_Alerts.hxx>
#include <BOPAlgo_BuilderSolid.hxx>
#include <BOPAlgo_CellsBuilder.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <IntCurvesFace_Intersector.hxx>
#include <IntTools_Context.hxx>
#include <Message_Report.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <Precision.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>

#include <BRepAlgoExt_Properties.hxx>

namespace BRepExt
{

namespace
{

// ---- The tolerances of the check, by name -------------------------------------------------
//
// T, "the partition's tolerance", is the largest of: OpenCascade's own Precision::Confusion()
// (1e-7), the fuzzy value the fuse ran at, and the largest tolerance of an edge or a face of
// the partition; but no more than tolerance::booleanCeiling. It is the width of the tube inside
// which two things may be taken for one, so nothing nearer than that to a boundary can say on which
// side it is. Vertex tolerances are left out of T: they are local (a ball round one point),
// and a probe at a vertex takes that vertex's own tolerance instead.

//! A point within kOnFactor x T of an operand's boundary is "on" it and decides nothing. Two,
//! because the point and the boundary may each be T off where they are meant to be.
constexpr double kOnFactor = 2;
//! A point inside a cell counts as inside only if it is at least kDepthFactor "on" distances
//! from the two faces it was found between; a cell thinner than that is "thin". Two: clear of
//! "on" with as much again to spare. (At four, a wall of 1e-6 between two cylinders was thin,
//! and a right Boolean refused.)
constexpr double kDepthFactor = 2;
// T is no more than tolerance::booleanCeiling, or than the largest tolerance an operand came
// with if that is more: ten times the document's tolerance of 1e-6. OpenCascade loosens an edge or
// a face as far as it needs to make two things one (a cylinder and its copy turned 1e-4
// rad: 1.18e-3; and as far as the caller's fuzzy value). Were T to follow, "on", "thin" and the
// allowance would grow with it and pass what is 1e-3 off. A section between curved faces in general
// position is loosened too (a torus and a ball: 6e-5) and is right, so a loose partition is
// not refused for being loose: it is checked at the ceiling, and what the fuse took for one
// must be one at that (check 4).

//! A face that the fuse took for a face of two operands is asked at up to kSharedSamples x
//! kSharedSamples places, or, where that grid misses it, along kSharedSamples lines across it
//! each way.
constexpr int kSharedSamples = 3;
//! How many faces of a cell are tried for a point inside it before the cell is called thin.
constexpr int kInsideTries = 6;
//! The fuzzy rung runs at kFuzzyFactor x the largest tolerance of the operands (not below
//! Precision::Confusion()): 1e-6 for operands at OpenCascade's default, which is the
//! document's tolerance. OpenCascade names the fuzzy value for "touching and near-coincident,
//! misaligned entities". Ten is one step up and no measurement: over 216 Booleans of the
//! reviewers' kinds this rung made none that the plain one had refused
//! It is never above tolerance::booleanFuzzyRelative of the
//! operands' diagonal, however loose they are.
constexpr double kFuzzyFactor = 10;
//! The cells of an operand must add up to the operand within
//!   kTube x T x (what bounds the faces that were cut)
//!     + tolerance::booleanQuadrature x (what those faces hold),
//! summed over the faces that had to be integrated at all (see Balance). The first term is
//! what a boundary that is T off can hold, on either side of it. The second is five times the
//! 1e-8 of a face that BRepAlgoExt_Properties.cxx integrates to; what was measured is far inside it
constexpr double kTube = 2;

using Hasher    = TopTools_ShapeMapHasher;
using ShapeList = NCollection_List<TopoDS_Shape>;

struct Clock
{
  static double now() { return booleanNow(); }
};

void appendOnce(ShapeSet& listed, ShapeList& list, const TopoDS_Shape& shape)
{
  if (listed.Contains(shape))
    return;
  listed.Add(shape);
  list.Append(shape);
}

double since(double start)
{
  return Clock::now() - start;
}

struct Freedom
{
  std::vector<std::pair<TopoDS_Shape, bool>> was;

  explicit Freedom(const BooleanRequest& request)
  {
    // Every part, not only the roots: a face or shell of an operand added to a new
    // container is frozen too.
    ShapeSet parts;
    for (const auto* group : {&request.sources, &request.tools})
      for (const auto& shape : *group)
        if (!shape.IsNull())
          TopExp::MapShapes(shape, parts);
    for (const auto& part : parts)
      was.emplace_back(part, part.TShape()->Free());
  }

  ~Freedom()
  {
    for (const auto& [shape, free] : was)
      shape.TShape()->Free(free);
  }
};

//! Why a rung does not give the result. Caught by the ladder, which goes on to the next.
struct Refused
{
  std::string why;
};

[[noreturn]] void refuse(const std::string& why)
{
  throw Refused{why};
}

const char* nameOf(Select select)
{
  switch (select)
  {
    case Select::Union:
      return "union";
    case Select::Difference:
      return "difference";
    case Select::Intersection:
    case Select::Common:
      return "intersection";
    case Select::Slice:
      return "slice";
    case Select::Imprint:
      return "imprint";
    case Select::Region:
      return "region";
    case Select::Behind:
      return "behind";
    case Select::Front:
      return "front";
    case Select::Auto:
      return "auto";
  }
  return "boolean";
}

//! How every refusal of a partition begins.
std::string inconsistent(const BooleanRequest& request)
{
  return std::string("the cells did not form a consistent ") + nameOf(request.select)
         + " partition: ";
}

std::string at(const gp_Pnt& p)
{
  char text[96];
  std::snprintf(text, sizeof text, "(%.4g, %.4g, %.4g)", p.X(), p.Y(), p.Z());
  return text;
}

//! 3 for solids, 2 for sheets, 1 for curves: the highest of what the operand holds. A body may
//! hold several kinds (a solid and a curve beside it); each is cut into cells of its own kind.
int dimensionOf(const char* operation, const TopoDS_Shape& shape)
{
  if (shape.IsNull())
    fail(operation, invalid_input, "an operand is empty");
  if (TopExp_Explorer(shape, TopAbs_SOLID).More())
    return 3;
  if (TopExp_Explorer(shape, TopAbs_FACE).More())
    return 2;
  if (TopExp_Explorer(shape, TopAbs_EDGE).More())
    return 1;
  fail(operation, invalid_input, "an operand has no solid, no face and no edge");
}

//! The solids of a shape, its faces that are of no solid, or its edges that are of no face.
TopExp_Explorer unitsOf(const TopoDS_Shape& shape, int dimension)
{
  return dimension == 3   ? TopExp_Explorer(shape, TopAbs_SOLID)
         : dimension == 2 ? TopExp_Explorer(shape, TopAbs_FACE, TopAbs_SOLID)
                          : TopExp_Explorer(shape, TopAbs_EDGE, TopAbs_FACE);
}

TopAbs_ShapeEnum unitType(int dimension)
{
  return dimension == 3 ? TopAbs_SOLID : dimension == 2 ? TopAbs_FACE : TopAbs_EDGE;
}

//! Whether two sets of shapes hold the same, in volume (area for sheets, length for curves),
//! without measuring what they have in common. Each side is added up face by face, one side
//! counted for and the other against; a face that both have, the same way round, and a face
//! between two cells of one side, cancel as they are counted, and are never integrated. What
//! is left is where the Boolean worked: the pieces of a face on one side and the face they
//! were cut from on the other. Those are integrated, all volumes about one origin: a flat
//! face bounded by lines and arcs by OpenCascade's plain Gauss rule, which is exact there and
//! quick; any other adaptively (BRepAlgoExt_Properties.cxx, to 1e-8 of the face), since with
//! splines in the surface or its outline the plain rule is off by more than a wrong result is
//! (2.4e-5 of a prism on a closed spline).
struct Balance
{
  gp_Pnt origin;
  NCollection_IndexedDataMap<TopoDS_Shape, int, Hasher>
    net[4]; // per dimension: face (or edge) -> count

  //! `twins`, if given, says of which face a face of `shape` is an exact copy: it is counted
  //! as that face, the way round the copy is.
  void add(const TopoDS_Shape&                                            shape,
           int                                                            dimension,
           int                                                            weight,
           const NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, Hasher>* twins = nullptr)
  {
    auto&      counts = net[dimension];
    const auto count  = [&](const TopoDS_Shape& found, int by) {
      const TopoDS_Shape* twin = twins ? twins->Seek(found) : nullptr;
      const TopoDS_Shape& part = twin ? *twin : found;
      if (int* seen = counts.ChangeSeek(part))
        *seen += by;
      else
        counts.Add(part, by);
    };
    if (dimension == 3)
    {
      // Of the solids alone: a face that is of no solid encloses nothing.
      for (TopExp_Explorer solid(shape, TopAbs_SOLID); solid.More(); solid.Next())
        for (TopExp_Explorer it(solid.Current(), TopAbs_FACE); it.More(); it.Next())
        {
          const auto way = it.Current().Orientation();
          if (way == TopAbs_FORWARD)
            count(it.Current(), weight);
          else if (way == TopAbs_REVERSED)
            count(it.Current(), -weight); // else it hangs inside and bounds nothing
        }
    }
    else if (dimension == 2)
    {
      for (TopExp_Explorer it(shape, TopAbs_FACE, TopAbs_SOLID); it.More(); it.Next())
        count(it.Current(), weight);
    }
    else
    {
      for (TopExp_Explorer it(shape, TopAbs_EDGE, TopAbs_FACE); it.More(); it.Next())
        count(it.Current(), weight);
    }
  }

  //! What the two sides differ by, and what that may be: see kTube and
  //! tolerance::booleanQuadrature.
  static bool plain(const TopoDS_Face& face)
  {
    if (BRepAdaptor_Surface(face, false).GetType() != GeomAbs_Plane)
      return false;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      const auto edge = TopoDS::Edge(it.Current());
      if (!BRep_Tool::Degenerated(edge) && BRepAdaptor_Curve(edge).GetType() > GeomAbs_Circle)
        return false;
    }
    return true;
  }

  void settle(int     dimension,
              double  tolerance,
              double& difference,
              double& allowed) const
  {
    difference = allowed = 0;
    const auto& counts   = net[dimension];
    for (int i = 1; i <= counts.Extent(); ++i)
    {
      if (counts(i) == 0)
        continue;
      GProp_GProps quick;
      double       held = 0;
      if (dimension == 1)
      {
        BRepGProp::LinearProperties(counts.FindKey(i), quick);
        held = quick.Mass();
        allowed += kTube * tolerance * 2; // its two ends
      }
      else
      {
        const auto face = TopoDS::Face(counts.FindKey(i).Oriented(TopAbs_FORWARD));
        // What bounds it, one measure down: its area for a volume, its outline for an area.
        // A plain face's area is what it holds, by the Gauss rule that is exact there; of
        // any other only the allowance and the integral's tolerance come from it.
        const bool isPlain = plain(face);
        double     bound;
        if (dimension == 3 && !isPlain)
        {
          bound = faceArea(face);
        }
        else
        {
          if (dimension == 3)
            BRepGProp::SurfaceProperties(face, quick);
          else
            BRepGProp::LinearProperties(face, quick);
          bound = quick.Mass();
        }
        allowed += kTube * tolerance * std::abs(bound);
        if (!isPlain)
        {
          const double area = bound;
          held = dimension == 3 ? faceMeasure(face, false, origin, nullptr, &area)
                                              : faceMeasure(face, dimension == 2, origin, nullptr);
        }
        else if (dimension == 3)
        {
          // On a plane (p - origin) . n is the same everywhere: its area times that.
          const auto   position = BRepAdaptor_Surface(face).Plane().Position();
          const gp_Vec normal =
            gp_Vec(position.XDirection()).Crossed(gp_Vec(position.YDirection()));
          held = quick.Mass() * gp_Vec(origin, position.Location()).Dot(normal) / 3;
        }
        else
        {
          GProp_GProps area;
          BRepGProp::SurfaceProperties(face, area);
          held = area.Mass();
        }
      }
      difference += counts(i) * held;
      allowed += tolerance::booleanQuadrature * std::abs(counts(i) * held);
    }
    difference = std::abs(difference);
  }
};

//! A volume, an area or a length by itself, where one is needed (a cell too thin to place).
double measureOf(const TopoDS_Shape& shape, int dimension, const gp_Pnt& origin)
{
  Balance alone;
  alone.origin = origin;
  alone.add(shape, dimension, 1);
  double total = 0;
  for (int i = 1; i <= alone.net[dimension].Extent(); ++i)
  {
    const auto& part = alone.net[dimension].FindKey(i);
    if (dimension == 1)
    {
      GProp_GProps quick;
      BRepGProp::LinearProperties(part, quick);
      total += quick.Mass();
    }
    else
      total +=
        alone.net[dimension](i)
        * faceMeasure(TopoDS::Face(part.Oriented(TopAbs_FORWARD)), dimension == 2, origin, nullptr);
  }
  return total;
}

//! The largest tolerance of an edge or a face.
double tubeOf(const TopoDS_Shape& shape)
{
  double   largest = 0;
  ShapeSet seen;
  TopExp::MapShapes(shape, TopAbs_EDGE, seen);
  for (int i = 1; i <= seen.Extent(); ++i)
    largest = std::max(largest, tolerance::own(seen(i)));
  seen.Clear();
  TopExp::MapShapes(shape, TopAbs_FACE, seen);
  for (int i = 1; i <= seen.Extent(); ++i)
    largest = std::max(largest, tolerance::own(seen(i)));
  return largest;
}

double vertexToleranceOf(const TopoDS_Shape& shape)
{
  double largest = 0;
  for (TopExp_Explorer it(shape, TopAbs_VERTEX); it.More(); it.Next())
    largest = std::max(largest, tolerance::own(it.Current()));
  return largest;
}

//! A face of a solid operand, with what a line is intersected with it by: made when the
//! first line comes near the face, and kept for the lines after it.
struct Skin
{
  TopoDS_Face                                face;
  Bnd_Box                                    box;
  std::unique_ptr<IntCurvesFace_Intersector> hits;
};

struct Lump
{
  TopoDS_Shape      solid;
  Bnd_Box           box;
  std::vector<Skin> skin;
};

//! The directions a point is looked out from, none along anything a model is likely to have.
const gp_Dir kOutlooks[] = {gp_Dir(0.5377, 0.2939, 0.7903),
                            gp_Dir(-0.6831, 0.6102, 0.4012),
                            gp_Dir(0.2319, -0.8547, 0.4645),
                            gp_Dir(-0.4129, -0.3391, -0.8453)};

//! Where a point is against one solid: a line from it crosses the solid's faces an odd number
//! of times if it is inside. A line that meets a face on its rim or only grazes it says
//! nothing, and another direction is taken, from `turn` on. Only faces whose boxes the line
//! goes through are asked; OpenCascade's own classifier prepares every face of the solid
//! first, which on a body of spline faces took longer than the Boolean.
TopAbs_State stateAgainst(Lump& lump, const gp_Pnt& point, double reach, int turn)
{
  for (int tried = 0; tried < 4; ++tried)
  {
    const gp_Lin line(point, kOutlooks[(turn + tried) % 4]);
    int          crossings = 0;
    bool         on = false, unclear = false;
    for (auto& skin : lump.skin)
    {
      if (skin.box.IsOut(line))
        continue;
      if (!skin.hits)
        skin.hits = std::make_unique<IntCurvesFace_Intersector>(skin.face, Precision::Confusion());
      skin.hits->Perform(line, -reach, RealLast());
      if (!skin.hits->IsDone())
      {
        unclear = true;
        break;
      }
      for (int i = 1; i <= skin.hits->NbPnt(); ++i)
      {
        const double along = skin.hits->WParameter(i);
        if (std::abs(along) <= reach)
          on = true;
        else if (along > 0)
        {
          unclear |=
            skin.hits->State(i) != TopAbs_IN || skin.hits->Transition(i) == IntCurveSurface_Tangent;
          ++crossings;
        }
      }
    }
    if (on)
      return TopAbs_ON;
    if (!unclear)
      return crossings % 2 ? TopAbs_IN : TopAbs_OUT;
  }
  return TopAbs_ON; // not known: it decides nothing
}

struct Piece
{
  TopoDS_Shape                shape;
  Bnd_Box                     box;
  std::shared_ptr<DistanceTo> distance; // made when first measured to
};

//! One operand, as the caller gave it: what is fused, and what everything is checked against.
struct Operand
{
  TopoDS_Shape original;
  TopoDS_Shape argument; // equivalent surface representation, otherwise the original
  NCollection_DataMap<TopoDS_Shape, ShapeList, Hasher> represented;
  bool                                                 tool      = false;
  int                                                  dimension = 0; // the highest it holds
  Bnd_Box                                              box;
  bool               has[4] = {}; // per dimension: it holds something of it
  std::vector<Lump>  lumps;       // its solids
  std::vector<Piece> pieces;      // its faces that are of no solid, its edges that are of no face
  ShapeSet           on; // everything of the partition that lies on this operand's boundary
};

struct Cell
{
  TopoDS_Shape      shape;
  int               dimension = 0;
  std::vector<char> owner; // per operand: the fuse made this cell out of it
  std::vector<char> in;    // per operand: a point inside the cell is inside (or on) it
  Bnd_Box           box;
  bool              thin = false;
};

struct Partition
{
  std::unique_ptr<BOPAlgo_CellsBuilder> builder;
  std::vector<Operand>                  operands;
  std::vector<Cell>                     cells;
  bool                                  represented = false;
  bool                                  journal     = false;
  double                                tolerance   = 0; // T
  Handle(IntTools_Context)              context     = new IntTools_Context;
  gp_Pnt origin; // about which every volume of the partition is taken
};

ShapeList imagesOf(BOPAlgo_CellsBuilder& builder, const TopoDS_Shape& shape)
{
  ShapeList images = builder.Modified(shape);
  if (images.IsEmpty() && !builder.IsDeleted(shape))
    images.Append(shape);
  return images;
}

// Compose a change of representation with the partition's history. Verification still
// measures and classifies against the original operands, including the shared-face check.
ShapeList imagesOf(Partition& partition, const TopoDS_Shape& original)
{
  if (!partition.represented)
    return imagesOf(*partition.builder, original);
  for (const auto& operand : partition.operands)
  {
    if (const auto* before = operand.represented.Seek(original))
    {
      ShapeList result;
      ShapeSet  seen;
      for (const auto& part : *before)
        for (const auto& image : imagesOf(*partition.builder, part))
          appendOnce(seen, result, image);
      return result;
    }
  }
  return imagesOf(*partition.builder, original);
}

// A complete sphere has one face; its pole/seam placement is arbitrary. Put the poles
// away from the near-coincident section, rather than accepting a lost cell at a pole.
bool wholeSphere(const TopoDS_Shape& shape, gp_Sphere& sphere, TopoDS_Face& face)
{
  if (shape.ShapeType() != TopAbs_SOLID)
    return false;
  TopExp_Explorer it(shape, TopAbs_FACE);
  if (!it.More())
    return false;
  face = TopoDS::Face(it.Current());
  it.Next();
  if (it.More())
    return false;
  BRepAdaptor_Surface surface(face);
  if (surface.GetType() != GeomAbs_Sphere
      || surface.LastUParameter() - surface.FirstUParameter()
           < 2 * M_PI - tolerance::booleanSurfaceParameter
      || surface.LastVParameter() - surface.FirstVParameter()
           < M_PI - tolerance::booleanSurfaceParameter)
    return false;
  sphere = surface.Sphere();
  return true;
}

bool spherePoles(std::vector<Operand>& operands)
{
  if (operands.size() != 2)
    return false;
  gp_Sphere   spheres[2];
  TopoDS_Face faces[2];
  if (!wholeSphere(operands[0].original, spheres[0], faces[0])
      || !wholeSphere(operands[1].original, spheres[1], faces[1]))
    return false;
  const gp_Vec between(spheres[0].Location(), spheres[1].Location());
  const double distance = between.Magnitude();
  if (distance <= Precision::Confusion()
      || distance
           > tolerance::booleanSphereRelative * std::min(spheres[0].Radius(), spheres[1].Radius()))
    return false;
  const gp_Vec along  = between.Normalized();
  const gp_Vec helper = std::abs(along.Z()) < 0.8 ? gp_Vec(0, 0, 1) : gp_Vec(0, 1, 0);
  const gp_Vec across = along.Crossed(helper).Normalized();
  const gp_Dir pole(along + across);
  for (size_t i = 0; i < 2; ++i)
  {
    auto& operand = operands[i];
    operand.argument =
      BRepPrimAPI_MakeSphere(gp_Ax2(spheres[i].Location(), pole, gp_Dir(along - across)),
                             spheres[i].Radius())
        .Shape();
    ShapeList solid, face;
    solid.Append(operand.argument);
    face.Append(TopExp_Explorer(operand.argument, TopAbs_FACE).Current());
    operand.represented.Bind(operand.original, solid);
    operand.represented.Bind(faces[i], face);
  }
  return true;
}

bool divideSeams(std::vector<Operand>& operands)
{
  bool changed = false;
  for (auto& operand : operands)
  {
    bool closed = false;
    for (TopExp_Explorer it(operand.original, TopAbs_FACE); it.More(); it.Next())
    {
      BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
      closed |= surface.GetType() == GeomAbs_Torus && surface.IsUClosed();
    }
    if (!closed)
      continue;
    BRepBuilderAPI_Copy            copier(operand.original, true, false);
    ShapeUpgrade_ShapeDivideClosed divide(copier.Shape());
    divide.SetNbSplitPoints(1);
    divide.Perform();
    const auto represented = divide.Result();
    if (represented.IsNull() || !BRepCheck_Analyzer(represented).IsValid())
      continue;
    operand.argument = represented;
    ShapeSet parts;
    TopExp::MapShapes(operand.original, parts);
    for (const auto& original : parts)
    {
      const auto now = divide.GetContext()->Apply(copier.ModifiedShape(original));
      ShapeList  images;
      if (now.ShapeType() == original.ShapeType())
        images.Append(now);
      else
        for (TopExp_Explorer it(now, original.ShapeType()); it.More(); it.Next())
          images.Append(it.Current());
      operand.represented.Bind(original, images);
    }
    changed = true;
  }
  return changed;
}

//! Warnings that are not refusals by themselves, because what they warn of is what the check
//! of the cells measures:
//! - a piece of a face that OpenCascade could not tell the way round of: the wrong way round
//!   it gives its cell a wrong volume;
//! - edges that no face was divided along: a line where two faces touch and do not cross
//!   divides nothing and is such an edge; one that should have divided a face leaves a cell
//!   that lies on both sides of an operand's boundary, which the probes find;
//! - faces that no solid was built from, where a sheet is among the operands: the part of a
//!   sheet that ends inside a solid bounds nothing. Between solids alone it stays a refusal.
bool checked(const Handle(Message_Alert)& alert, bool sheets)
{
  return alert->IsKind(STANDARD_TYPE(BOPAlgo_AlertUnableToOrientTheShape))
         || alert->IsKind(STANDARD_TYPE(BOPAlgo_AlertFaceBuilderUnusedEdges))
         || (sheets && alert->IsKind(STANDARD_TYPE(BOPAlgo_AlertSolidBuilderUnusedFaces)));
}

//! OpenCascade's alerts of one gravity, by name; of its warnings those that `checked` leaves.
std::string alerts(const BOPAlgo_CellsBuilder& builder,
                   Message_Gravity             gravity,
                   bool                        sheets = false)
{
  std::string names;
  for (const auto& alert : builder.GetReport()->GetAlerts(gravity))
  {
    if (gravity == Message_Warning && checked(alert, sheets))
      continue;
    const std::string name = alert->DynamicType()->Name();
    if (names.find(name) != std::string::npos)
      continue;
    if (!names.empty())
      names += ", ";
    names += name;
  }
  return names;
}

//! Where `point` is against the original operand: inside, outside, or within `reach` of its
//! boundary. A sheet or a curve has no inside: a point is on it or off it.
//! For a point inside a solid cell only the operand's solids are asked (`solidsOnly`).
TopAbs_State stateAgainst(Operand&      operand,
                          const gp_Pnt& point,
                          double        reach,
                          bool          solidsOnly = false,
                          int           turn       = 0)
{
  Bnd_Box here;
  here.Add(point);
  here.Enlarge(reach);
  if (operand.box.IsOut(here))
    return TopAbs_OUT;
  try
  {
    TopAbs_State best = TopAbs_OUT;
    for (auto& lump : operand.lumps)
    {
      if (lump.box.IsOut(here))
        continue;
      const auto state = stateAgainst(lump, point, reach, turn);
      if (state == TopAbs_IN)
        return TopAbs_IN;
      if (state != TopAbs_OUT)
        best = TopAbs_ON;
    }
    if (best == TopAbs_ON || solidsOnly || operand.pieces.empty())
      return best;
    const auto vertex = BRepBuilderAPI_MakeVertex(point).Vertex();
    for (auto& piece : operand.pieces)
    {
      if (piece.box.IsOut(here))
        continue;
      if (!piece.distance)
        piece.distance = std::make_shared<DistanceTo>(piece.shape);
      const auto apart = piece.distance->from(vertex);
      if (apart && *apart <= reach)
        return TopAbs_ON;
    }
    return TopAbs_OUT;
  }
  catch (const Standard_Failure&)
  {
    return TopAbs_ON; // not known: it decides nothing
  }
}

//! A point of a face, an edge or a vertex, away from what bounds it. False if none was found.
bool pointOn(const TopoDS_Shape& element, const Handle(IntTools_Context)& context, gp_Pnt& point)
{
  try
  {
    if (element.ShapeType() == TopAbs_VERTEX)
    {
      point = BRep_Tool::Pnt(TopoDS::Vertex(element));
      return true;
    }
    if (element.ShapeType() == TopAbs_EDGE)
    {
      const auto edge = TopoDS::Edge(element);
      if (BRep_Tool::Degenerated(edge))
        return false;
      BRepAdaptor_Curve curve(edge);
      // Off the middle, where symmetric shapes have their edges and seams.
      point = curve.Value(curve.FirstParameter()
                          + (curve.LastParameter() - curve.FirstParameter()) * 0.3819660112501051);
      return true;
    }
    gp_Pnt2d uv;
    return BOPTools_AlgoTools3D::PointInFace(TopoDS::Face(element), point, uv, context) == 0;
  }
  catch (const Standard_Failure&)
  {
    return false;
  }
}

//! A point inside a solid cell and how deep it is: from a point of one of its faces, straight
//! in, half way to the first thing the cell has there. The count of what the line meets must
//! be odd, or the line was not going in (a cell thinner than `floor`, or a grazing line).
bool pointInside(const TopoDS_Shape&             cell,
                 double                          floor,
                 const Handle(IntTools_Context)& context,
                 gp_Pnt&                         point,
                 double&                         depth)
{
  depth = 0;
  try
  {
    // The faces of the cell, each asked only by a line through its box, and given its
    // intersector then (as an operand's skin is): the cell's other faces cannot be met.
    std::vector<Skin> faces;
    for (TopExp_Explorer it(cell, TopAbs_FACE); it.More(); it.Next())
    {
      faces.emplace_back();
      faces.back().face = TopoDS::Face(it.Current());
      BRepBndLib::Add(it.Current(), faces.back().box, false);
    }
    int tried = 0;
    for (TopExp_Explorer it(cell, TopAbs_FACE); it.More() && tried < kInsideTries; it.Next())
    {
      const auto face = TopoDS::Face(it.Current());
      if (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED)
        continue;
      gp_Pnt   from;
      gp_Pnt2d uv;
      if (BOPTools_AlgoTools3D::PointInFace(face, from, uv, context) != 0)
        continue;
      BRepAdaptor_Surface surface(face);
      gp_Pnt              same;
      gp_Vec              du, dv;
      surface.D1(uv.X(), uv.Y(), same, du, dv);
      gp_Vec outward = du.Crossed(dv);
      if (outward.Magnitude() < tolerance::fineAngle)
        continue;
      if (face.Orientation() == TopAbs_REVERSED)
        outward.Reverse();
      ++tried;
      const gp_Lin line(from, gp_Dir(outward.Reversed()));
      int          crossings = 0;
      double       nearest   = RealLast();
      for (auto& skin : faces)
      {
        if (skin.box.IsOut(line))
          continue;
        if (!skin.hits)
          skin.hits =
            std::make_unique<IntCurvesFace_Intersector>(skin.face, Precision::Confusion());
        skin.hits->Perform(line, floor, RealLast());
        if (!skin.hits->IsDone())
          continue; // as IntCurvesFace_ShapeIntersector counts it: no points
        crossings += skin.hits->NbPnt();
        for (int i = 1; i <= skin.hits->NbPnt(); ++i)
          nearest = std::min(nearest, skin.hits->WParameter(i));
      }
      if (crossings % 2 == 0)
        continue;
      if (nearest / 2 <= depth)
        continue;
      depth = nearest / 2;
      point = from.Translated(gp_Vec(gp_Dir(outward.Reversed())) * depth);
      if (depth > 1000 * floor)
        break; // deep enough; no need to look for deeper
    }
  }
  catch (const Standard_Failure&)
  {
  }
  return depth > 0;
}

// ---- The fuse ------------------------------------------------------------------------------

Operand operandOf(const BooleanRequest& request, const TopoDS_Shape& shape, bool tool)
{
  Operand operand;
  operand.original  = shape;
  operand.tool      = tool;
  operand.dimension = dimensionOf(request.operation, shape);
  BRepBndLib::Add(shape, operand.box, false);
  for (auto it = unitsOf(shape, 3); it.More(); it.Next())
  {
    operand.lumps.emplace_back();
    auto& lump = operand.lumps.back();
    lump.solid = it.Current();
    BRepBndLib::Add(it.Current(), lump.box, false);
    for (TopExp_Explorer face(it.Current(), TopAbs_FACE); face.More(); face.Next())
    {
      const auto way = face.Current().Orientation();
      if (way != TopAbs_FORWARD && way != TopAbs_REVERSED)
        continue; // hangs inside, bounds nothing
      lump.skin.emplace_back();
      lump.skin.back().face = TopoDS::Face(face.Current());
      BRepBndLib::Add(face.Current(), lump.skin.back().box, false);
    }
  }
  for (const int dimension : {2, 1})
    for (auto it = unitsOf(shape, dimension); it.More(); it.Next())
    {
      operand.pieces.emplace_back();
      operand.pieces.back().shape = it.Current();
      BRepBndLib::Add(it.Current(), operand.pieces.back().box, false);
    }
  return operand;
}

//! Where two walls that were swept along the same direction cross, they cross in straight
//! lines along it, from one rim to the other. OpenCascade looks for such a line by walking
//! over both surfaces and loses it where the walls cross at a shallow angle (a prism and a
//! copy moved a little); where a rim meets the other wall it finds without trouble. So the
//! lines are made here, from where the rims of each wall meet the other, and given to the
//! fuse beside the operands: it cuts both walls along them.
ShapeList sweptCrossings(const std::vector<Operand>& operands, double crossingTolerance = 1e-4)
{
  struct Wall
  {
    TopoDS_Face face;
    gp_Dir      along;
    Bnd_Box     box;
  };

  std::vector<std::vector<Wall>> walls(operands.size());
  for (size_t o = 0; o < operands.size(); ++o)
    for (TopExp_Explorer it(operands[o].original, TopAbs_FACE); it.More(); it.Next())
    {
      BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
      if (surface.GetType() != GeomAbs_SurfaceOfExtrusion)
        continue;
      walls[o].push_back({TopoDS::Face(it.Current()), surface.Direction(), {}});
      BRepBndLib::Add(it.Current(), walls[o].back().box, false);
    }
  ShapeList lines;
  // A rim that lies in the other wall's surface does not cross it: the two walls are one
  // there (two prisms on one outline), which OpenCascade finds by itself. Over the wall's own
  // parameters: its surface runs on without end along the sweep, and a search spread over all
  // of that took milliseconds a point; a rim beyond the wall does not meet it either.
  const auto lies = [](const BRepAdaptor_Curve& curve, GeomAPI_ProjectPointOnSurf& foot) {
    for (const double at : {0.21, 0.5, 0.83})
    {
      foot.Perform(curve.Value(curve.FirstParameter()
                               + at * (curve.LastParameter() - curve.FirstParameter())));
      if (foot.NbPoints() == 0 || foot.LowerDistance() > Precision::Confusion())
        return false;
    }
    return true;
  };
  const auto footOn = [](const Wall& wall, GeomAPI_ProjectPointOnSurf& foot) {
    double u0, u1, v0, v1;
    BRepTools::UVBounds(wall.face, u0, u1, v0, v1);
    foot.Init(BRep_Tool::Surface(wall.face), u0, u1, v0, v1, Precision::Confusion());
  };
  // A line lies along both walls, each of which it crosses on a line of its own along the
  // sweep: its curve on each is that line between the feet of its ends, given with it, so that
  // the fuse does not project it. Where an end or its middle is not on the wall the fuse does.
  const auto curveOn = [](const TopoDS_Edge&          line,
                          const gp_Pnt&               low,
                          const gp_Pnt&               high,
                          const Wall&                 wall,
                          GeomAPI_ProjectPointOnSurf& foot) {
    const double tolerance = BRep_Tool::Tolerance(line);
    gp_Pnt2d     ends[2];
    int          k = 0;
    for (const gp_Pnt& end : {low, high})
    {
      foot.Perform(end);
      if (foot.NbPoints() == 0 || foot.LowerDistance() > tolerance)
        return;
      double u, v;
      foot.LowerDistanceParameters(u, v);
      ends[k++] = gp_Pnt2d(u, v);
    }
    double first, last;
    BRep_Tool::Range(line, first, last);
    NCollection_Array1<gp_Pnt2d> poles(1, 2);
    poles(1) = ends[0];
    poles(2) = ends[1];
    NCollection_Array1<double> knots(1, 2);
    knots(1) = first;
    knots(2) = last;
    NCollection_Array1<int> multiplicities(1, 2);
    multiplicities.Init(2);
    const Handle(Geom2d_BSplineCurve) onWall =
      new Geom2d_BSplineCurve(poles, knots, multiplicities, 1);
    const gp_Pnt2d middle = onWall->Value((first + last) / 2);
    if (BRep_Tool::Surface(wall.face)
          ->Value(middle.X(), middle.Y())
          .Distance(gp_Pnt((low.XYZ() + high.XYZ()) / 2))
        > tolerance)
      return;
    BRep_Builder().UpdateEdge(line, onWall, wall.face, tolerance);
  };
  const auto meet = [&lies](const Wall&                 rims,
                            const Wall&                 wall,
                            GeomAPI_ProjectPointOnSurf& foot,
                            std::vector<gp_Pnt>&        points) {
    IntCurvesFace_Intersector hits(wall.face, Precision::Confusion());
    for (TopExp_Explorer it(rims.face, TopAbs_EDGE); it.More(); it.Next())
    {
      const auto edge = TopoDS::Edge(it.Current());
      if (BRep_Tool::Degenerated(edge))
        continue;
      Handle(BRepAdaptor_Curve) curve = new BRepAdaptor_Curve(edge);
      if (lies(*curve, foot))
        continue;
      hits.Perform(curve, curve->FirstParameter(), curve->LastParameter());
      if (!hits.IsDone())
        continue;
      for (int i = 1; i <= hits.NbPnt(); ++i)
        points.push_back(hits.Pnt(i));
    }
  };
  for (size_t a = 0; a < operands.size(); ++a)
    for (size_t b = a + 1; b < operands.size(); ++b)
      for (const auto& one : walls[a])
        for (const auto& other : walls[b])
        {
          if (one.box.IsOut(other.box)
              || !one.along.IsParallel(other.along, BRepExt::tolerance::angle))
            continue;
          GeomAPI_ProjectPointOnSurf onOne, onOther;
          footOn(one, onOne);
          footOn(other, onOther);
          std::vector<gp_Pnt> points;
          meet(one, other, onOther, points);
          meet(other, one, onOne, points);
          // Points over each other along the sweep are one line.
          const gp_Vec      along(one.along);
          std::vector<char> used(points.size(), 0);
          for (size_t i = 0; i < points.size(); ++i)
          {
            if (used[i])
              continue;
            gp_Pnt low = points[i], high = points[i];
            for (size_t j = i + 1; j < points.size(); ++j)
            {
              const gp_Vec between(points[i], points[j]);
              if (used[j] || between.Crossed(along).Magnitude() > crossingTolerance)
                continue;
              used[j] = 1;
              if (gp_Vec(low, points[j]).Dot(along) < 0)
                low = points[j];
              if (gp_Vec(high, points[j]).Dot(along) > 0)
                high = points[j];
            }
            if (low.Distance(high) <= BRepExt::tolerance::sewing)
              continue;
            const TopoDS_Edge line = BRepBuilderAPI_MakeEdge(low, high).Edge();
            curveOn(line, low, high, one, onOne);
            curveOn(line, low, high, other, onOther);
            lines.Append(line);
          }
        }
  return lines;
}

//! One way of running the fuse.
struct Step
{
  const char* name;
  double      fuzzy; // OpenCascade's fuzzy value
};

//! `lines` are given to the fuse beside the operands: they are cut with them and cut them, and
//! are none of them, so they make no cell.
void fuse(Partition&            partition,
          const BooleanRequest& request,
          const Step&           step,
          const ShapeList&      lines,
          std::vector<Operand>* prepared = nullptr)
{
  partition.journal  = request.journal;
  const double fuzzy = step.fuzzy;
  const double setup = booleanTimingEnabled() ? booleanNow() : 0;
  if (prepared)
    partition.operands = std::move(*prepared);
  else
  {
    for (const auto& shape : request.sources)
      partition.operands.push_back(operandOf(request, shape, false));
    for (const auto& shape : request.tools)
      partition.operands.push_back(operandOf(request, shape, true));
  }
  // Boundary images belong to this rung, while boxes/classifiers belong to immutable inputs.
  for (auto& operand : partition.operands)
  {
    operand.on.Clear();
    partition.represented |= !operand.represented.IsEmpty();
  }
  auto&   operands = partition.operands;
  Bnd_Box all;
  for (const auto& operand : operands)
    all.Add(operand.box);
  if (!all.IsVoid() && !all.IsOpen())
    partition.origin = gp_Pnt((all.CornerMin().XYZ() + all.CornerMax().XYZ()) / 2);

  if (booleanTimingEnabled())
    booleanPhase("operandOf", since(setup));
  BooleanPhase metadata("fuseMetadata");
  partition.builder = std::make_unique<BOPAlgo_CellsBuilder>();
  auto&     builder = *partition.builder;
  ShapeList arguments;
  for (const auto& operand : operands)
    arguments.Append(operand.argument.IsNull() ? operand.original : operand.argument);
  for (const auto& line : lines)
    arguments.Append(line);
  builder.SetArguments(arguments);
  // The operands are the document's, and stay as they were.
  builder.SetNonDestructive(true);
  builder.SetRunParallel(false);
  builder.SetUseOBB(true);
  builder.SetGlue(static_cast<BOPAlgo_GlueEnum>(request.glue));
  if (fuzzy > 0)
    builder.SetFuzzyValue(fuzzy);
  {
    BooleanPhase phase("perform");
    const auto   indicator = progress(request.cancellable);
    builder.Perform(indicator.IsNull() ? Message_ProgressRange() : indicator->Start());
    checkCancelled(request.operation, request.cancellable);
  }
  if (builder.HasErrors())
    refuse(inconsistent(request) + "OpenCascade could not intersect the operands ("
           + alerts(builder, Message_Fail) + ")");
  // A warning is OpenCascade saying that it left something out or patched something up:
  // an intersection it gave up on, edges a face was not built from, faces a solid was not.
  bool sheets = false;
  for (const auto& operand : operands)
    sheets |= operand.dimension < 3;
  if (const auto warnings = alerts(builder, Message_Warning, sheets); !warnings.empty())
    refuse(inconsistent(request) + "OpenCascade intersected the operands with reservations ("
           + warnings + ")");
  builder.AddAllToResult(0, false);

  partition.tolerance = std::max({Precision::Confusion(), fuzzy, tubeOf(builder.GetAllParts())});
  double ceiling      = tolerance::booleanCeiling;
  for (const auto& operand : operands)
    ceiling = std::max(ceiling, tubeOf(operand.original));
  partition.tolerance = std::min(partition.tolerance, ceiling);

  ShapeSet index;
  for (size_t o = 0; o < operands.size(); ++o)
  {
    checkCancelled(request.operation, request.cancellable);
    auto& operand = operands[o];
    for (const int dimension : {3, 2, 1})
      for (auto it = unitsOf(operand.original, dimension); it.More(); it.Next())
      {
        operand.has[dimension] = true;
        for (const auto& image : imagesOf(partition, it.Current()))
        {
          if (image.ShapeType() != unitType(dimension))
            continue;
          const int found = index.Add(image);
          if (found > static_cast<int>(partition.cells.size()))
          {
            partition.cells.emplace_back();
            auto& cell     = partition.cells.back();
            cell.shape     = image;
            cell.dimension = dimension;
            cell.owner.assign(operands.size(), 0);
            cell.in.assign(operands.size(), 0);
          }
          partition.cells[found - 1].owner[o] = 1;
        }
      }
    // What of the partition lies on this operand's boundary: its faces as they are now
    // (for a curve its edges), with everything that bounds them.
    const auto mark = [&](const TopoDS_Shape& piece) {
      for (const auto& image : imagesOf(partition, piece))
      {
        TopExp::MapShapes(image, TopAbs_FACE, operand.on);
        TopExp::MapShapes(image, TopAbs_EDGE, operand.on);
        TopExp::MapShapes(image, TopAbs_VERTEX, operand.on);
      }
    };
    for (TopExp_Explorer it(operand.original, TopAbs_FACE); it.More(); it.Next())
      mark(it.Current());
    for (auto it = unitsOf(operand.original, 1); it.More(); it.Next())
      mark(it.Current());
  }
}

// ---- The check -----------------------------------------------------------------------------
//
// 1. Every cell is in or out of each operand: one point inside the cell is classified against
//    the ORIGINAL operands, and where the fuse has its own word on that (which operands it
//    made the cell from), the two must agree.
// 2. The cells of each operand add up to the operand, in volume (area for a sheet, length for
//    a curve). Summed over two operands this is V(A) + V(B) = V(A or B) + V(A and B), and
//    V(A - B) = V(A) - V(A and B) is then how the difference is chosen.
// 3. No cell says otherwise on its boundary: a point of each face, edge and vertex of a cell
//    that is not itself on an operand's boundary must be on the side of that operand that the
//    cell is on. This is what finds operands that overlap and were not cut by each other
//    (where 1 and 2 hold: each cell is one whole operand), since a boundary that crosses a
//    cell unnoticed has that cell's faces, edges or corners on both of its sides.

void verify(Partition& partition, const BooleanRequest& request, Rung& rung)
{
  auto&             operands     = partition.operands;
  const double      reach        = kOnFactor * partition.tolerance;
  const std::string inconsistent = BRepExt::inconsistent(request);

  BooleanPhase check("verify");
  // 1.
  const double                       interiorStart = booleanTimingEnabled() ? booleanNow() : 0;
  std::vector<std::array<double, 4>> unverified(operands.size(), std::array<double, 4>{});
  for (auto& cell : partition.cells)
  {
    checkCancelled(request.operation, request.cancellable);
    // A cell that is a solid (a face, an edge) of an operand as it was is inside that
    // operand, and outside every operand whose box it does not reach: nothing to ask.
    std::vector<char> settled(operands.size(), 0);
    BRepBndLib::Add(cell.shape, cell.box, false);
    cell.box.Enlarge(reach);
    const Bnd_Box& room = cell.box;
    bool           open = false;
    for (size_t o = 0; o < operands.size(); ++o)
    {
      checkCancelled(request.operation, request.cancellable);
      if (operands[o].dimension < cell.dimension)
        settled[o] = 1;
      else if (operands[o].box.IsOut(room))
        settled[o] = 1;
      else if (cell.owner[o])
      {
        for (auto it = unitsOf(operands[o].original, cell.dimension); it.More() && !settled[o];
             it.Next())
          if (it.Current().IsSame(cell.shape))
            settled[o] = 2;
      }
      open |= !settled[o];
    }
    if (!open)
    {
      for (size_t o = 0; o < operands.size(); ++o)
        cell.in[o] = settled[o] == 2;
      for (size_t o = 0; o < operands.size(); ++o)
        if (operands[o].has[cell.dimension] && cell.in[o] != cell.owner[o])
          refuse(inconsistent + "a cell of an operand lies wholly outside the box of that operand");
      continue;
    }
    gp_Pnt     inside;
    double     depth = RealLast();
    auto       timed = Clock::now();
    const bool found =
      cell.dimension == 3
        ? pointInside(cell.shape, kDepthFactor * reach, partition.context, inside, depth)
        : pointOn(cell.shape, partition.context, inside);
    rung.insideMs += since(timed);
    timed     = Clock::now();
    cell.thin = !found || depth < kDepthFactor * reach;
    std::vector<TopAbs_State> states(operands.size(), TopAbs_OUT);
    if (!cell.thin)
    {
      for (size_t o = 0; o < operands.size(); ++o)
      {
        checkCancelled(request.operation, request.cancellable);
        if (settled[o])
        {
          states[o] = settled[o] == 2 ? TopAbs_IN : TopAbs_OUT;
          continue;
        }
        states[o] = stateAgainst(operands[o], inside, reach, cell.dimension == 3);
        ++rung.interior;
        // A sheet lies in a solid also where it lies in its skin; a solid's own cell
        // has to be clear of every boundary to say anything.
        if (cell.dimension == 3 && states[o] == TopAbs_ON)
          cell.thin = true;
      }
    }
    rung.classifyMs += since(timed);
    if (cell.thin)
    {
      // Too thin to say where it is: the fuse's word stands, and what such cells hold
      // is counted against what the identities allow.
      cell.in           = cell.owner;
      const double held = std::abs(measureOf(cell.shape, cell.dimension, partition.origin));
      for (size_t o = 0; o < operands.size(); ++o)
        if (cell.owner[o])
          unverified[o][cell.dimension] += held;
      continue;
    }
    for (size_t o = 0; o < operands.size(); ++o)
    {
      checkCancelled(request.operation, request.cancellable);
      cell.in[o] = operands[o].dimension >= cell.dimension && states[o] != TopAbs_OUT;
      // The fuse has a word on it only for operands that hold something of the cell's kind.
      if (!operands[o].has[cell.dimension] || cell.in[o] == cell.owner[o])
        continue;
      // Looked at once more, from other directions, before it is a refusal.
      const auto again = stateAgainst(operands[o], inside, reach, cell.dimension == 3, 2);
      if (again == TopAbs_ON || (again != TopAbs_OUT) == (cell.owner[o] != 0))
      {
        // Undecided, a solid cell is taken on the fuse's word like one too thin to place.
        if (again == TopAbs_ON && cell.dimension == 3)
          unverified[o][3] += std::abs(measureOf(cell.shape, 3, partition.origin));
        cell.in[o] = cell.owner[o];
        continue;
      }
      refuse(inconsistent + "the cell at " + at(inside) + " is "
             + (cell.in[o] ? "inside" : "outside") + " an operand that the fuse "
             + (cell.owner[o] ? "made it from" : "did not make it from"));
    }
  }

  if (booleanTimingEnabled())
    booleanPhase("interior", since(interiorStart));
  const double balanceStart = booleanTimingEnabled() ? booleanNow() : 0;
  // 2.
  for (size_t o = 0; o < operands.size(); ++o)
    for (const int dimension : {3, 2, 1})
    {
      if (!operands[o].has[dimension])
        continue;
      Balance balance;
      balance.origin = partition.origin;
      for (const auto& cell : partition.cells)
        if (cell.owner[o] && cell.dimension == dimension)
          balance.add(cell.shape, dimension, 1);
      balance.add(operands[o].original, dimension, -1);
      double     residual, slack;
      const auto timed = Clock::now();
      balance.settle(dimension, partition.tolerance, residual, slack);
      rung.measureMs += since(timed);
      if (residual / std::max(slack, 1e-300) > rung.residual / std::max(rung.slack, 1e-300))
      {
        rung.residual = residual;
        rung.slack    = slack;
      }
      char numbers[160];
      if (residual > slack)
      {
        std::snprintf(numbers,
                      sizeof numbers,
                      "the cells of an operand hold %.9g more or less than the operand, where %.3g "
                      "is allowed (missing or excess material)",
                      residual,
                      slack);
        refuse(inconsistent + numbers);
      }
      if (unverified[o][dimension] > slack)
      {
        std::snprintf(numbers,
                      sizeof numbers,
                      "cells thinner than the tolerance hold %.9g of an operand",
                      unverified[o][dimension]);
        refuse(inconsistent + numbers);
      }
    }

  if (booleanTimingEnabled())
    booleanPhase("balance", since(balanceStart));
  const double sharedStart = booleanTimingEnabled() ? booleanNow() : 0;
  // 4. What the fuse took for one is one. A face (an edge) of the partition that the fuse
  //    has as the image of a face (an edge) of two operands lies on both of those, within
  //    the same "on" distance, wherever it is asked. OpenCascade takes two things for one
  //    within whatever tolerance it has loosened them to, or the caller's fuzzy value: two
  //    cylinders of 5000 moved 1e-4 across each other came back with one rim for the two
  //    lower ones and no bottom between them.
  for (const bool faces : {true, false})
  {
    const auto type = faces ? TopAbs_FACE : TopAbs_EDGE;
    NCollection_IndexedDataMap<TopoDS_Shape, std::vector<std::pair<size_t, TopoDS_Shape>>, Hasher>
      origins;
    for (size_t o = 0; o < operands.size(); ++o)
      for (TopExp_Explorer it(operands[o].original, type); it.More(); it.Next())
        for (const auto& image : imagesOf(partition, it.Current()))
        {
          if (image.ShapeType() != type)
            continue;
          if (!origins.Contains(image))
            origins.Add(image, {});
          origins.ChangeFromKey(image).emplace_back(o, it.Current());
        }
    for (int index = 1; index <= origins.Extent(); ++index)
    {
      const auto& from    = origins(index);
      bool        several = false;
      for (const auto& origin : from)
        several |= origin.first != from.front().first;
      if (!several)
        continue;
      std::vector<gp_Pnt> points;
      if (faces)
      {
        const auto            face = TopoDS::Face(origins.FindKey(index));
        BRepAdaptor_Surface   surface(face);
        const double          u0 = surface.FirstUParameter(), u1 = surface.LastUParameter();
        const double          v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
        std::vector<gp_Pnt2d> inside;
        if (!Precision::IsInfinite(u0) && !Precision::IsInfinite(u1) && !Precision::IsInfinite(v0)
            && !Precision::IsInfinite(v1))
          for (int i = 0; i < kSharedSamples; ++i)
            for (int j = 0; j < kSharedSamples; ++j)
            {
              const gp_Pnt2d uv(u0 + (u1 - u0) * (i + 0.5) / kSharedSamples,
                                v0 + (v1 - v0) * (j + 0.5) / kSharedSamples);
              if (partition.context->IsPointInFace(face, uv))
                inside.push_back(uv);
            }
        // A face the grid misses (a ring, a slot) is asked along lines across it instead,
        // and one with no point to ask at all is not taken for one on trust.
        if (inside.empty())
          inside = parametersAcross(face, kSharedSamples, partition.context);
        if (inside.empty())
          refuse(inconsistent
                 + "two faces were taken for one that has no point inside it to check");
        for (const auto& uv : inside)
          points.push_back(surface.Value(uv.X(), uv.Y()));
      }
      else
      {
        const auto edge = TopoDS::Edge(origins.FindKey(index));
        if (BRep_Tool::Degenerated(edge))
          continue;
        BRepAdaptor_Curve curve(edge);
        for (const double along : {0.21, 0.5, 0.83})
          points.push_back(curve.Value(curve.FirstParameter()
                                       + along * (curve.LastParameter() - curve.FirstParameter())));
      }
      double inputTolerance = Precision::Confusion();
      for (const auto& ancestor : from)
        inputTolerance = std::max(inputTolerance, tubeOf(ancestor.second));
      const double sharedReach = kOnFactor * inputTolerance;
      for (const auto& point : points)
        for (const auto& [o, original] : from)
        {
          if (original.IsSame(origins.FindKey(index)))
            continue; // it is that face (edge) itself
          double apart = RealLast();
          if (faces)
          {
            GeomAPI_ProjectPointOnSurf foot(point, BRep_Tool::Surface(TopoDS::Face(original)));
            if (foot.NbPoints() > 0)
              apart = foot.LowerDistance();
          }
          else
          {
            BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(point).Vertex(),
                                                original);
            if (distance.IsDone())
              apart = distance.Value();
          }
          // A section may acquire a loose tolerance without making two distinct
          // input boundaries the same. A shared image must agree at the INPUTS'
          // tolerance, rather than at the tolerance the fuse raised to merge them.
          if (apart <= sharedReach)
            continue;
          char numbers[220];
          std::snprintf(numbers,
                        sizeof numbers,
                        "two %s were taken for one that are %.3g apart at %s, where the input "
                        "boundaries allow %.3g",
                        faces ? "faces" : "edges",
                        apart,
                        at(point).c_str(),
                        sharedReach);
          refuse(inconsistent + numbers);
        }
    }
  }

  if (booleanTimingEnabled())
    booleanPhase("shared", since(sharedStart));
  BooleanPhase boundary("boundary");
  // 3.
  ShapeSet                              elements;
  std::vector<std::vector<signed char>> known; // per element and operand: -1 not asked yet
  std::vector<gp_Pnt>                   points;
  std::vector<signed char>              located; // -1 not asked yet, 0 no point, 1 a point
  std::vector<Bnd_Box>                  reaches; // the box of each element
  const auto                            stateOf = [&](int element, size_t o) {
    auto& cached = known[element - 1][o];
    if (cached >= 0)
      return static_cast<TopAbs_State>(cached);
    const auto&  shape = elements(element);
    TopAbs_State state = TopAbs_ON;
    if (!operands[o].on.Contains(shape))
    {
      auto& box = reaches[element - 1];
      if (box.IsVoid())
        BRepBndLib::Add(shape, box, false);
      // Nothing of what the operand's box does not reach is inside the operand.
      if (operands[o].box.IsOut(box))
      {
        cached = static_cast<signed char>(TopAbs_OUT);
        return TopAbs_OUT;
      }
      if (located[element - 1] < 0)
        located[element - 1] = pointOn(shape, partition.context, points[element - 1]) ? 1 : 0;
      if (located[element - 1] == 1)
      {
        const double own = shape.ShapeType() == TopAbs_VERTEX ? tolerance::own(shape) : 0;
        const auto   timed = Clock::now();
        state =
          stateAgainst(operands[o], points[element - 1], std::max(reach, kOnFactor * own), false);
        rung.classifyMs += since(timed);
        ++rung.probes;
      }
    }
    cached = static_cast<signed char>(state);
    return state;
  };
  for (const auto& cell : partition.cells)
  {
    std::vector<size_t> candidates;
    for (size_t o = 0; o < operands.size(); ++o)
      if (operands[o].dimension >= cell.dimension
          && (cell.in[o] || !operands[o].box.IsOut(cell.box)))
        candidates.push_back(o);
    for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
    {
      if ((type == TopAbs_FACE && cell.dimension < 3)
          || (type == TopAbs_EDGE && cell.dimension < 2))
        continue;
      for (TopExp_Explorer it(cell.shape, type); it.More(); it.Next())
      {
        const int element = elements.Add(it.Current());
        if (element > static_cast<int>(known.size()))
        {
          known.emplace_back(operands.size(), -1);
          points.emplace_back();
          located.push_back(-1);
          reaches.emplace_back();
        }
        for (const size_t o : candidates)
        {
          const auto state = stateOf(element, o);
          const bool wrong = cell.in[o] ? state == TopAbs_OUT : state == TopAbs_IN;
          if (!wrong)
            continue;
          // Looked at once more, from other directions, before it is a refusal.
          if (located[element - 1] != 1)
          {
            // Outside the operand's box there is nothing to look at again.
            pointOn(elements(element), partition.context, points[element - 1]);
          }
          else if (stateAgainst(operands[o], points[element - 1], reach, false, 2) != state)
          {
            known[element - 1][o] = static_cast<signed char>(TopAbs_ON);
            continue;
          }
          refuse(inconsistent + "a cell that is " + (cell.in[o] ? "inside" : "outside")
                 + " an operand has a "
                 + (type == TopAbs_FACE   ? "face"
                    : type == TopAbs_EDGE ? "edge"
                                          : "corner")
                 + " at " + at(points[element - 1]) + " that is "
                 + (cell.in[o] ? "outside" : "inside")
                 + " it: the operands overlap there and were not cut by each other");
        }
      }
    }
  }
}

// ---- Choosing and putting together -----------------------------------------------------------

//! Takes out of a body's faces the edges of intersection that lie inside a face and divide nothing:
//! where two operands touch along a line in a face, the partition keeps that line in the face as a
//! wire of INTERNAL edges. It has the face on both sides, so the face's points are the same without
//! it, and it is an edge on one face in the shell. In place: these wires are the partition's own,
//! and an edge an operand brought with it is left alone.
void dropIntersectionWithin(const TopoDS_Shape& body, const BooleanRequest& request)
{
  NCollection_IndexedMap<TopoDS_Shape, Hasher> given;
  for (const auto* group : {&request.sources, &request.tools})
    for (const auto& operand : *group)
      TopExp::MapShapes(operand, TopAbs_EDGE, given);
  BRep_Builder                                 builder;
  NCollection_IndexedMap<TopoDS_Shape, Hasher> faces;
  TopExp::MapShapes(body, TopAbs_FACE, faces);
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    std::vector<TopoDS_Shape> within;
    for (TopoDS_Iterator wire(faces(i)); wire.More(); wire.Next())
    {
      if (wire.Value().ShapeType() != TopAbs_WIRE)
        continue;
      bool internal = false, hanging = true;
      for (TopoDS_Iterator edge(wire.Value()); edge.More() && hanging; edge.Next())
      {
        internal = true;
        hanging  = edge.Value().Orientation() == TopAbs_INTERNAL && !given.Contains(edge.Value());
      }
      if (internal && hanging)
        within.push_back(wire.Value());
    }
    if (within.empty())
      continue;
    TopoDS_Shape face = faces(i);
    const bool   free = face.Free();
    face.Free(true);
    for (const auto& wire : within)
      builder.Remove(face, wire);
    face.Free(free);
  }
}

//! A solid without the faces and shells that hang inside it and bound nothing.
TopoDS_Shape withoutInternals(const TopoDS_Shape& solid)
{
  bool any = false;
  for (TopoDS_Iterator shell(solid); shell.More() && !any; shell.Next())
  {
    any =
      shell.Value().Orientation() == TopAbs_INTERNAL || shell.Value().ShapeType() != TopAbs_SHELL;
    if (!any)
      for (TopoDS_Iterator face(shell.Value()); face.More() && !any; face.Next())
        any = face.Value().Orientation() == TopAbs_INTERNAL;
  }
  if (!any)
    return solid;
  BRep_Builder builder;
  TopoDS_Solid clean;
  builder.MakeSolid(clean);
  for (TopoDS_Iterator shell(solid); shell.More(); shell.Next())
  {
    if (shell.Value().Orientation() == TopAbs_INTERNAL || shell.Value().ShapeType() != TopAbs_SHELL)
      continue;
    TopoDS_Shell kept;
    builder.MakeShell(kept);
    bool any = false;
    for (TopoDS_Iterator face(shell.Value()); face.More(); face.Next())
      if (face.Value().Orientation() != TopAbs_INTERNAL)
      {
        builder.Add(kept, face.Value());
        any = true;
      }
    if (!any)
      continue; // a shell of hanging faces alone
    kept.Closed(true);
    kept.Orientation(shell.Value().Orientation());
    builder.Add(clean, kept);
  }
  return clean;
}

TopoDS_Shape solidFrom(const ShapeList& faces, const BooleanRequest& request)
{
  BOPAlgo_BuilderSolid maker;
  maker.SetShapes(faces);
  const auto indicator = progress(request.cancellable);
  maker.Perform(indicator.IsNull() ? Message_ProgressRange() : indicator->Start());
  checkCancelled(request.operation, request.cancellable);
  if (maker.HasErrors() || maker.HasWarnings() || maker.Areas().Extent() != 1)
    refuse("the chosen cells did not close into one solid");
  return maker.Areas().First();
}

//! Cells that touch across a face become one solid, without that face.
std::vector<TopoDS_Shape> joined(const std::vector<TopoDS_Shape>& cells,
                                 const BooleanRequest&            request)
{
  std::vector<TopoDS_Shape> solids;
  if (cells.empty())
    return solids;
  BRep_Builder    builder;
  TopoDS_Compound all;
  builder.MakeCompound(all);
  for (const auto& cell : cells)
    builder.Add(all, cell);
  ShapeList blocks;
  BOPTools_AlgoTools::MakeConnexityBlocks(all, TopAbs_FACE, TopAbs_SOLID, blocks);
  for (const auto& block : blocks)
  {
    checkCancelled(request.operation, request.cancellable);
    NCollection_IndexedDataMap<TopoDS_Shape, int, Hasher> uses;
    int                                                   count = 0;
    TopoDS_Shape                                          only;
    for (TopoDS_Iterator solid(block); solid.More(); solid.Next(), ++count)
    {
      only = solid.Value();
      for (TopoDS_Iterator shell(solid.Value()); shell.More(); shell.Next())
      {
        if (shell.Value().Orientation() == TopAbs_INTERNAL
            || shell.Value().ShapeType() != TopAbs_SHELL)
          continue;
        for (TopoDS_Iterator face(shell.Value()); face.More(); face.Next())
        {
          if (face.Value().Orientation() == TopAbs_INTERNAL)
            continue;
          if (int* seen = uses.ChangeSeek(face.Value()))
            ++*seen;
          else
            uses.Add(face.Value(), 1);
        }
      }
    }
    if (count == 1)
    {
      solids.push_back(withoutInternals(only));
      continue;
    }
    ShapeList outer;
    for (int i = 1; i <= uses.Extent(); ++i)
      if (uses(i) == 1)
        outer.Append(uses.FindKey(i));
    solids.push_back(solidFrom(outer, request));
  }
  return solids;
}

//! Faces as sheets: those in `groups[i]` make one shell (a face alone stays a face).
TopoDS_Shape sheetOf(const std::vector<TopoDS_Shape>& faces)
{
  if (faces.size() == 1)
    return faces.front();
  BRep_Builder builder;
  TopoDS_Shell shell;
  builder.MakeShell(shell);
  for (const auto& face : faces)
    builder.Add(shell, face);
  BOPTools_AlgoTools::OrientFacesOnShell(shell);
  return shell;
}

//! Faces (or edges) gathered into what hangs together, across every shared edge (or end) that
//! is not in `walls`.
std::vector<std::vector<TopoDS_Shape>> connected(
  const std::vector<TopoDS_Shape>&                parts,
  TopAbs_ShapeEnum                                across,
  const std::function<bool(const TopoDS_Shape&)>& wall)
{
  std::vector<int> group(parts.size());
  std::iota(group.begin(), group.end(), 0);
  const std::function<int(int)> root = [&](int i) {
    return group[i] == i ? i : group[i] = root(group[i]);
  };
  NCollection_DataMap<TopoDS_Shape, int, Hasher> first;
  for (size_t i = 0; i < parts.size(); ++i)
    for (TopExp_Explorer it(parts[i], across); it.More(); it.Next())
    {
      if (wall(it.Current()))
        continue;
      if (const int* other = first.Seek(it.Current()))
        group[root(static_cast<int>(i))] = root(*other);
      else
        first.Bind(it.Current(), static_cast<int>(i));
    }
  std::vector<std::vector<TopoDS_Shape>> groups;
  std::vector<int>                       slot(parts.size(), -1);
  for (size_t i = 0; i < parts.size(); ++i)
  {
    const int r = root(static_cast<int>(i));
    if (slot[r] < 0)
    {
      slot[r] = static_cast<int>(groups.size());
      groups.emplace_back();
    }
    groups[slot[r]].push_back(parts[i]);
  }
  return groups;
}

TopoDS_Shape wireOf(const std::vector<TopoDS_Shape>& edges)
{
  BRep_Builder builder;
  TopoDS_Wire  wire;
  builder.MakeWire(wire);
  for (const auto& edge : edges)
    builder.Add(wire, edge);
  return wire;
}

//! A solid with a sheet, joined: what the sheet closes off together with the solid's skin is
//! added to the solid (the enclosed regions of a solid and a sheet body). The pockets are built
//! from the solid's faces seen from outside and the sheet's pieces outside it, seen from both
//! sides.
std::vector<TopoDS_Shape> enclosed(Partition&            partition,
                                   const BooleanRequest& request,
                                   Rung&                 rung,
                                   bool                  regions = false)
{
  ShapeList                 faces;
  std::vector<TopoDS_Shape> solids;
  for (const auto& cell : partition.cells)
  {
    if (cell.dimension == 3)
    {
      solids.push_back(cell.shape);
      continue;
    }
    bool inside = false;
    for (size_t o = 0; o < partition.operands.size(); ++o)
      inside |= partition.operands[o].dimension == 3 && cell.in[o];
    if (inside)
      continue;
    faces.Append(cell.shape.Oriented(TopAbs_FORWARD));
    faces.Append(cell.shape.Oriented(TopAbs_REVERSED));
  }
  // The skin of the solids: the faces that one cell alone has.
  for (const auto& whole : joined(solids, request))
    for (TopExp_Explorer it(whole, TopAbs_FACE); it.More(); it.Next())
      faces.Append(it.Current().Reversed());
  BOPAlgo_BuilderSolid maker;
  maker.SetShapes(faces);
  const auto indicator = progress(request.cancellable);
  maker.Perform(indicator.IsNull() ? Message_ProgressRange() : indicator->Start());
  checkCancelled(request.operation, request.cancellable);
  if (maker.HasErrors())
    refuse("the sheet and the solid's skin could not be built into regions");
  const double reach   = kOnFactor * partition.tolerance;
  int          pockets = 0;
  for (const auto& pocket : maker.Areas())
  {
    gp_Pnt inside;
    double depth = 0;
    if (!pointInside(pocket, kDepthFactor * reach, partition.context, inside, depth))
      continue;
    bool outside = true;
    for (auto& operand : partition.operands)
      if (operand.dimension == 3)
      {
        outside &= stateAgainst(operand, inside, reach) == TopAbs_OUT;
        ++rung.interior;
      }
    if (!outside || measureOf(pocket, 3, partition.origin) <= 0)
      continue;
    solids.push_back(pocket);
    ++pockets;
  }
  if (pockets == 0 && (!regions || solids.empty()))
    refuse(
      regions
        ? "the sheets enclose no closed region"
        : "the sheet closes off no region together with the solid, so there is nothing to join");
  return regions ? solids : joined(solids, request);
}

// ---- Sides of a sheet ------------------------------------------------------------------------
//
// A sheet tool given a side of material (face normals pointing away from the material, or
// into it) bounds what it cuts through: a cell of a source lies
// behind it or in front of it. The side is read where the cell meets the sheet, against the
// normal of the tool's own face there, as the caller gave it.

enum class Side
{
  None,
  Behind,
  Front,
  Both
};

Side joinSides(Side a, Side b)
{
  if (a == Side::None)
    return b;
  if (b == Side::None || a == b)
    return a;
  return Side::Both;
}

//! The faces of the partition that lie on the sheet tools, each with the tool's face it is a
//! piece of; and the edges of those pieces, each with every tool face it borders.
struct SheetFaces
{
  NCollection_DataMap<TopoDS_Shape, TopoDS_Face, Hasher>              faceOf;
  NCollection_DataMap<TopoDS_Shape, std::vector<TopoDS_Face>, Hasher> edgeOf;

  explicit SheetFaces(Partition& partition)
  {
    for (const auto& operand : partition.operands)
    {
      if (!operand.tool)
        continue;
      for (TopExp_Explorer face(operand.original, TopAbs_FACE); face.More(); face.Next())
        for (const auto& image : imagesOf(partition, face.Current()))
        {
          if (!faceOf.IsBound(image))
            faceOf.Bind(image, TopoDS::Face(face.Current()));
          for (TopExp_Explorer edge(image, TopAbs_EDGE); edge.More(); edge.Next())
          {
            if (!edgeOf.IsBound(edge.Current()))
              edgeOf.Bind(edge.Current(), {});
            auto& faces = edgeOf.ChangeFind(edge.Current());
            if (std::none_of(faces.begin(), faces.end(), [&](const TopoDS_Face& known) {
                  return known.IsSame(face.Current());
                }))
              faces.push_back(TopoDS::Face(face.Current()));
          }
        }
    }
  }
};

//! The normal of the tool's face `face` (as the tool has it) at its point nearest `point`.
bool toolNormal(const TopoDS_Face&              face,
                const gp_Pnt&                   point,
                const Handle(IntTools_Context)& context,
                gp_Dir&                         normal)
{
  GeomAPI_ProjectPointOnSurf& projector = context->ProjPS(face);
  projector.Perform(point);
  if (!projector.IsDone() || projector.NbPoints() == 0)
    return false;
  double u = 0, v = 0;
  projector.LowerDistanceParameters(u, v);
  BRepAdaptor_Surface surface(face, false);
  gp_Pnt              at;
  gp_Vec              du, dv;
  surface.D1(u, v, at, du, dv);
  gp_Vec n = du.Crossed(dv);
  if (n.Magnitude() < tolerance::fineAngle)
    return false;
  if (face.Orientation() == TopAbs_REVERSED)
    n.Reverse();
  normal = gp_Dir(n);
  return true;
}

//! Which side of the sheet tools `cell` lies on. A solid cell is bounded by pieces of the
//! sheet: it lies on the side its outward normal there points away from. A face of a sheet
//! source meets the sheet along edges: it lies on the side it leaves the sheet towards, which
//! a point of it just off such an edge says when it is clear of "on" (kOnFactor x T).
Side sideOf(const Cell& cell, const SheetFaces& sheets, Partition& partition)
{
  Side side = Side::None;
  if (cell.dimension == 3)
  {
    for (TopExp_Explorer it(cell.shape, TopAbs_FACE); it.More(); it.Next())
    {
      const auto* tool = sheets.faceOf.Seek(it.Current());
      if (!tool)
        continue;
      const auto face = TopoDS::Face(it.Current());
      // A piece of the sheet that hangs inside the cell divides nothing.
      if (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED)
        return Side::Both;
      gp_Pnt   point;
      gp_Pnt2d uv;
      if (BOPTools_AlgoTools3D::PointInFace(face, point, uv, partition.context) != 0)
        continue;
      BRepAdaptor_Surface surface(face);
      gp_Pnt              at;
      gp_Vec              du, dv;
      surface.D1(uv.X(), uv.Y(), at, du, dv);
      gp_Vec outward = du.Crossed(dv);
      gp_Dir normal;
      if (outward.Magnitude() < tolerance::fineAngle
          || !toolNormal(*tool, point, partition.context, normal))
        continue;
      if (face.Orientation() == TopAbs_REVERSED)
        outward.Reverse();
      side = joinSides(side, outward.Dot(gp_Vec(normal)) > 0 ? Side::Behind : Side::Front);
    }
    return side;
  }
  const double reach = kOnFactor * partition.tolerance;
  const auto   face  = TopoDS::Face(cell.shape.Oriented(TopAbs_FORWARD));
  for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
  {
    const auto* tools = sheets.edgeOf.Seek(it.Current());
    const auto  edge  = TopoDS::Edge(it.Current());
    if (!tools || BRep_Tool::Degenerated(edge) || BRep_Tool::IsClosed(edge, face))
      continue;
    double first = 0, last = 0;
    BRep_Tool::Range(edge, first, last);
    const double t = first + (last - first) * 0.3819660112501051; // off the middle, as pointOn
    gp_Pnt2d     uv;
    gp_Pnt       off;
    try
    {
      if (BOPTools_AlgoTools3D::PointNearEdge(edge, face, t, uv, off, partition.context) != 0)
        continue;
    }
    catch (const Standard_Failure&)
    {
      continue;
    }
    const gp_Pnt on   = BRepAdaptor_Curve(edge).Value(t);
    Side         here = Side::None;
    for (const auto& tool : *tools)
    {
      gp_Dir normal;
      if (!toolNormal(tool, on, partition.context, normal))
        continue;
      const double height = gp_Vec(on, off).Dot(gp_Vec(normal));
      if (std::abs(height) > reach)
        here = joinSides(here, height > 0 ? Side::Front : Side::Behind);
    }
    // Where the sheet folds along the edge its two faces may disagree: that edge says nothing.
    if (here != Side::Both)
      side = joinSides(side, here);
  }
  return side;
}

struct Chosen
{
  BooleanResult::Decided    decided = BooleanResult::Decided::Asked;
  std::vector<TopoDS_Shape> bodies; // solids, sheets or wires
  std::vector<TopoDS_Shape> cells;  // what they were put together from
  bool                      merged = false;
};

Chosen choose(Partition& partition, const BooleanRequest& request, Rung& rung)
{
  const auto& operands = partition.operands;
  if (request.select == Select::Region)
  {
    for (const auto& operand : operands)
      if (operand.dimension < 2)
        fail(request.operation,
             unsupported_geometry,
             "regions require solids or sheets, not curves");
    Chosen regions;
    if (operands.front().dimension == 3)
    {
      for (const auto& cell : partition.cells)
        if (cell.dimension == 3)
          regions.bodies.push_back(cell.shape);
    }
    else
      regions.bodies = enclosed(partition, request, rung, true);
    regions.cells = regions.bodies;
    return regions;
  }
  int sourceDimension = 0, toolDimension = 0, lowest = 3;
  for (const auto& operand : operands)
  {
    int& dimension = operand.tool ? toolDimension : sourceDimension;
    if (dimension != 0 && dimension != operand.dimension)
      fail(request.operation,
           invalid_input,
           "the bodies on one side of a Boolean must all be solids, or all sheets, or all curves");
    dimension = operand.dimension;
    lowest    = std::min(lowest, operand.dimension);
  }
  Select select = request.select;
  Chosen chosen;
  if (select == Select::Auto)
  {
    if (sourceDimension != 3 || toolDimension != 3)
      fail(request.operation, invalid_input, "an automatic Boolean takes solids on both sides");
    // The cells are verified: each knows which operands it lies in, so the clash is read
    // off them instead of being asked of the kernel again.
    bool     overlap = false;
    ShapeSet sourceFaces;
    for (const auto& cell : partition.cells)
    {
      if (cell.dimension != 3)
        continue;
      bool inSource = false, inTool = false;
      for (size_t o = 0; o < operands.size(); ++o)
        (operands[o].tool ? inTool : inSource) |= cell.in[o] != 0;
      overlap |= inSource && inTool;
      if (inSource && !inTool)
        for (TopExp_Explorer it(cell.shape, TopAbs_FACE); it.More(); it.Next())
          sourceFaces.Add(it.Current());
    }
    bool touching = false;
    for (const auto& cell : partition.cells)
    {
      if (overlap || touching)
        break;
      if (cell.dimension != 3)
        continue;
      bool inSource = false, inTool = false;
      for (size_t o = 0; o < operands.size(); ++o)
        (operands[o].tool ? inTool : inSource) |= cell.in[o] != 0;
      if (!inTool || inSource)
        continue;
      for (TopExp_Explorer it(cell.shape, TopAbs_FACE); it.More() && !touching; it.Next())
        touching = sourceFaces.Contains(it.Current());
    }
    if (!overlap && !touching)
    {
      chosen.decided = BooleanResult::Decided::None;
      return chosen;
    }
    select         = overlap ? Select::Difference : Select::Union;
    chosen.decided = overlap ? BooleanResult::Decided::Difference : BooleanResult::Decided::Union;
  }
  if (select == Select::Behind || select == Select::Front)
  {
    if (toolDimension != 2 || sourceDimension < 2)
      fail(request.operation,
           invalid_input,
           "a side of a sheet is chosen with sheets as tools, from solids or sheets");
    const SheetFaces          sheets(partition);
    const Side                wanted = select == Select::Behind ? Side::Behind : Side::Front;
    std::vector<TopoDS_Shape> kept;
    for (const auto& cell : partition.cells)
    {
      checkCancelled(request.operation, request.cancellable);
      if (cell.dimension != sourceDimension)
        continue;
      bool source = false, onTool = false;
      for (size_t o = 0; o < operands.size(); ++o)
        (operands[o].tool ? onTool : source) |=
          (operands[o].tool ? cell.in[o] : cell.owner[o]) != 0;
      if (!source)
        continue;
      if (onTool)
        fail(request.operation,
             unsupported_geometry,
             "a piece of the body lies on the sheet itself, which has no side for it");
      const Side side = sideOf(cell, sheets, partition);
      if (side == Side::None)
        fail(request.operation,
             no_solution,
             "the sheet does not divide the body: a piece of it does not meet the sheet, so it is "
             "on neither side");
      if (side == Side::Both)
        fail(request.operation,
             no_solution,
             "a piece of the body meets the sheet from both sides: the sheet does not divide it, "
             "or its faces are not turned alike");
      if (side != wanted)
        continue;
      kept.push_back(cell.shape);
      chosen.cells.push_back(cell.shape);
    }
    if (sourceDimension == 3)
    {
      chosen.bodies = joined(kept, request);
      chosen.merged = !kept.empty();
    }
    else
    {
      const auto never = [](const TopoDS_Shape&) { return false; };
      for (const auto& group : connected(kept, TopAbs_EDGE, never))
        chosen.bodies.push_back(sheetOf(group));
    }
    return chosen;
  }

  // A sheet has no inside to take away: it divides (a solid sectioned by a sheet).
  if (select == Select::Difference && toolDimension != 0 && toolDimension < sourceDimension)
    select = Select::Slice;

  if (select == Select::Union && toolDimension != 0 && toolDimension != sourceDimension)
  {
    if (sourceDimension != 3 || toolDimension != 2)
      fail(
        request.operation,
        unsupported_geometry,
        "only a sheet can be joined to a solid, not the other way round, and no curve to either");
    chosen.bodies = enclosed(partition, request, rung);
    chosen.merged = true;
    chosen.cells  = chosen.bodies; // the pockets have no cells of the fuse to be held against
    return chosen;
  }
  // What operands of different kinds have in common is of the lowest kind among them.
  const bool lowestOnly = select == Select::Intersection || select == Select::Common;

  std::vector<TopoDS_Shape>              cells[4];
  std::vector<std::vector<TopoDS_Shape>> perSource(operands.size());
  for (const auto& cell : partition.cells)
  {
    if (lowestOnly && cell.dimension != lowest)
      continue;
    bool source = false, inSource = false, inTool = false, inAll = true;
    for (size_t o = 0; o < operands.size(); ++o)
    {
      checkCancelled(request.operation, request.cancellable);
      if (!operands[o].tool)
      {
        source |= cell.owner[o] != 0;
        inSource |= cell.in[o] != 0;
      }
      else
        inTool |= cell.in[o] != 0;
      inAll &= cell.in[o] != 0;
    }
    bool take = false;
    switch (select)
    {
      case Select::Region:
      case Select::Behind:
      case Select::Front:
      case Select::Auto:
        break; // handled before cell selection
      case Select::Union:
        take = true;
        break;
      case Select::Difference:
        take = source && !inTool;
        break;
      case Select::Intersection:
        take = inSource && inTool;
        break;
      case Select::Common:
        take = inAll;
        break;
      case Select::Slice:
      case Select::Imprint:
        take = source;
        break;
    }
    if (!take)
      continue;
    cells[cell.dimension].push_back(cell.shape);
    chosen.cells.push_back(cell.shape);
    // A cell that two sources share is of both, and each comes back whole.
    if (select == Select::Imprint && cell.dimension == 3)
    {
      bool first = true;
      for (size_t o = 0; o < operands.size(); ++o)
        if (!operands[o].tool && cell.owner[o])
        {
          perSource[o].push_back(cell.shape);
          if (!first)
            chosen.cells.push_back(cell.shape);
          first = false;
        }
    }
  }

  // Where a sheet (a curve) is divided: along every edge (at every end) that a tool has or
  // that the sources did not have before, which is one that the fuse cut into them.
  ShapeSet had;
  for (const auto& operand : operands)
  {
    if (operand.tool)
      continue;
    for (const auto type : {TopAbs_EDGE, TopAbs_VERTEX})
      for (TopExp_Explorer it(operand.original, type); it.More(); it.Next())
        for (const auto& image : imagesOf(partition, it.Current()))
          had.Add(image);
  }
  const auto wall = [&](const TopoDS_Shape& edge) {
    if (!had.Contains(edge))
      return true;
    for (const auto& operand : operands)
      if (operand.tool && operand.on.Contains(edge))
        return true;
    return false;
  };
  const auto never = [](const TopoDS_Shape&) { return false; };
  if (select == Select::Slice)
  {
    for (const auto& cell : cells[3])
      chosen.bodies.push_back(withoutInternals(cell));
  }
  else if (select == Select::Imprint)
  {
    for (const auto& own : perSource)
      for (const auto& solid : joined(own, request))
        chosen.bodies.push_back(solid);
    chosen.merged = true;
  }
  else
  {
    chosen.bodies = joined(cells[3], request);
    chosen.merged = !cells[3].empty();
  }
  const std::function<bool(const TopoDS_Shape&)> divides =
    select == Select::Slice ? std::function<bool(const TopoDS_Shape&)>(wall) : never;
  for (const auto& group : connected(cells[2], TopAbs_EDGE, divides))
    chosen.bodies.push_back(sheetOf(group));
  for (const auto& group : connected(cells[1], TopAbs_VERTEX, divides))
    chosen.bodies.push_back(wireOf(group));
  return chosen;
}

//! The faces of the operands as the bodies have them now.
using Families = NCollection_DataMap<TopoDS_Shape, ShapeList, Hasher>;

Families familiesOf(Partition&                                                     partition,
                    const TopoDS_Shape&                                            shape,
                    const NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, Hasher>& roots)
{
  ShapeSet present;
  TopExp::MapShapes(shape, present);
  Families families;
  for (const auto& operand : partition.operands)
    for (auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
      for (TopExp_Explorer it(operand.original, type); it.More(); it.Next())
      {
        if (type != TopAbs_FACE && !partition.journal)
          continue;
        if (families.IsBound(it.Current()))
          continue;
        ShapeList now;
        ShapeSet  listed;
        for (const auto& image : imagesOf(partition, it.Current()))
        {
          const auto* owned = roots.Seek(image);
          const auto& final = owned ? *owned : image;
          if (final.ShapeType() == it.Current().ShapeType() && present.Contains(final))
            appendOnce(listed, now, final);
        }
        families.Bind(it.Current(), now);
      }
  return families;
}

//! Pieces of one original face that lie side by side in one body of the result are that face
//! again, and an edge that ran through to a new vertex is one edge again. Faces of different
//! origin keep the edge between them, and every vertex that an operand had stays.
//!
//! It is done on a copy of the body, geometry and all. The body shares every face that the
//! fuse left alone with an operand, and ShapeUpgrade_UnifySameDomain changes what it works on
//! in place where a surface is closed (it moved the curves of a torus operand's face on its
//! surface, which left the operand invalid and with it the body that was returned). `twins`
//! receives, for each face of the copy, the face it is a copy of.
TopoDS_Shape reunite(Partition&                                               partition,
                     const TopoDS_Shape&                                      shape,
                     const std::vector<TopoDS_Shape>&                         bodies,
                     Families&                                                families,
                     NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, Hasher>& twins,
                     Handle(BRepTools_History)&                               transformation,
                     std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>>&      displayCopies)
{
  NCollection_DataMap<TopoDS_Shape, int, Hasher> bodyOf;
  for (size_t body = 0; body < bodies.size(); ++body)
    for (TopExp_Explorer it(bodies[body], TopAbs_FACE); it.More(); it.Next())
      bodyOf.Bind(it.Current(),
                  bodyOf.IsBound(it.Current()) ? -1 : static_cast<int>(body)); // -1: of two bodies
  NCollection_DataMap<TopoDS_Shape, std::vector<int>, Hasher> origins;
  int                                                         origin = 0;
  for (Families::Iterator it(families); it.More(); it.Next(), ++origin)
    for (const auto& piece : it.Value())
    {
      if (piece.ShapeType() != TopAbs_FACE)
        continue;
      if (!origins.IsBound(piece))
        origins.Bind(piece, {});
      origins.ChangeFind(piece).push_back(origin);
    }
  Adjacency beside;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, beside);
  std::vector<TopoDS_Shape> kept;
  bool                      any = false;
  for (int e = 1; e <= beside.Extent(); ++e)
  {
    const auto& faces   = beside(e);
    bool        release = false;
    const int*  in      = faces.IsEmpty() ? nullptr : bodyOf.Seek(faces.First());
    const int*  other   = faces.IsEmpty() ? nullptr : bodyOf.Seek(faces.Last());
    if (faces.Extent() == 2 && !faces.First().IsSame(faces.Last()) && in && other && *in >= 0
        && *in == *other)
    {
      const auto* a = origins.Seek(faces.First());
      const auto* b = origins.Seek(faces.Last());
      if (a && b)
        for (const int one : *a)
          release |= std::find(b->begin(), b->end(), one) != b->end();
    }
    if (release)
      any = true;
    else
      kept.push_back(beside.FindKey(e));
  }
  if (!any)
    return shape;

  BRepBuilderAPI_Copy copier(shape, true, true);
  ShapeSet            present;
  TopExp::MapShapes(shape, TopAbs_VERTEX, present);
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    if (!twins.IsBound(copier.ModifiedShape(it.Current())))
    {
      const auto* original = twins.Seek(it.Current());
      twins.Bind(copier.ModifiedShape(it.Current()), original ? *original : it.Current());
    }
  ShapeUpgrade_UnifySameDomain unify(copier.Shape(), true, true, true);
  for (const auto& edge : kept)
    unify.KeepShape(copier.ModifiedShape(edge));
  for (const auto& operand : partition.operands)
    for (TopExp_Explorer it(operand.original, TopAbs_VERTEX); it.More(); it.Next())
      for (const auto& image : imagesOf(partition, it.Current()))
        if (present.Contains(image))
          unify.KeepShape(copier.ModifiedShape(image));
  unify.Build();
  const auto& united = unify.Shape();
  if (partition.journal)
  {
    ShapeSet original, final;
    TopExp::MapShapes(shape, original);
    TopExp::MapShapes(united, final);
    for (int i = 1; i <= original.Extent(); ++i)
    {
      const auto& before = original(i);
      if (!BRepTools_History::IsSupportedType(before))
        continue;
      const auto  copy     = copier.ModifiedShape(before);
      const auto& modified = unify.History()->Modified(copy);
      if (modified.IsEmpty() && !unify.History()->IsRemoved(copy) && final.Contains(copy))
      {
        transformation->AddModified(before, copy);
        displayCopies.emplace_back(before, copy);
      }
      for (const auto& after : modified)
        if (final.Contains(after))
          transformation->AddModified(before, after);
      for (const auto& after : unify.History()->Generated(copy))
        if (final.Contains(after))
          transformation->AddGenerated(before, after);
    }
  }

  for (Families::Iterator it(families); it.More(); it.Next())
  {
    ShapeList now;
    ShapeSet  listed;
    for (const auto& original : it.Value())
    {
      const TopoDS_Shape piece   = copier.ModifiedShape(original);
      const auto&        changed = unify.History()->Modified(piece);
      if (changed.IsEmpty() && !unify.History()->IsRemoved(piece))
        appendOnce(listed, now, piece);
      for (const auto& face : changed)
        if (face.ShapeType() == it.Key().ShapeType())
          appendOnce(listed, now, face);
    }
    it.ChangeValue() = now;
  }
  return united;
}

BooleanResult build(Partition& partition, const BooleanRequest& request, Rung& rung)
{
  BooleanPhase buildPhase("build");
  auto         chosen = choose(partition, request, rung);
  if (chosen.decided == BooleanResult::Decided::None)
  {
    BooleanResult none;
    none.decided = chosen.decided;
    return none;
  }
  for (const auto& body : chosen.bodies)
    if (body.ShapeType() == TopAbs_SOLID)
      dropIntersectionWithin(body, request);
  BRep_Builder    builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, Hasher> roots, twins;
  for (auto& body : chosen.bodies)
  {
    // Add freezes its child. Give the result its own root TShape before adding it:
    // restoring an input's Free flag must never unfreeze a child of this result.
    auto owned = body.EmptyCopied();
    owned.Free(true);
    for (TopoDS_Iterator child(body, false, false); child.More(); child.Next())
      builder.Add(owned, child.Value());
    builder.Add(compound, owned);
    // A one-face sheet is its root. Account for the owned face in history and
    // in the balance, just as for faces copied during same-domain reunion.
    if (body.ShapeType() == TopAbs_FACE)
    {
      roots.Bind(body, owned);
      twins.Bind(owned, body);
    }
    body = owned;
  }

  BooleanResult result;
  result.shape   = compound;
  result.decided = chosen.decided;
  auto families  = familiesOf(partition, compound, roots);
  // What was put together is the cells that were chosen, and is a body.
  const auto sound = [&](const TopoDS_Shape& shape) {
    if (!accepts(shape))
      return false;
    for (const int dimension : {3, 2})
    {
      Balance balance;
      balance.origin = partition.origin;
      balance.add(shape, dimension, 1, &twins);
      for (const auto& cell : chosen.cells)
        balance.add(cell, dimension, -1);
      double difference, allowed;
      balance.settle(dimension, partition.tolerance, difference, allowed);
      if (difference > allowed)
        return false;
    }
    return true;
  };
  if (chosen.merged ? !sound(compound) : !accepts(compound))
    refuse(inconsistent(request) + "the chosen cells did not make a valid body of their volume");
  Handle(BRepTools_History) transformation = new BRepTools_History;
  if (request.select != Select::Imprint && request.select != Select::Region)
  {
    // Tidying is not owed: where the pieces of a face do not go back into one face that
    // is sound, they stay pieces, and the body is right as it is.
    auto                                               reunited = families;
    std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> displayCopies;
    const auto                                         tidy =
      reunite(partition, compound, chosen.bodies, reunited, twins, transformation, displayCopies);
    if (!tidy.IsSame(compound) && sound(tidy))
    {
      result.shape         = tidy;
      result.displayCopies = std::move(displayCopies);
      families             = reunited;
    }
  }
  for (Families::Iterator it(families); it.More(); it.Next())
  {
    const auto& now = it.Value();
    if (now.Extent() == 1 && now.First().IsSame(it.Key()))
      continue; // in the result as it was
    if (it.Key().ShapeType() == TopAbs_FACE)
      result.faces.Bind(it.Key(), now);
  }
  if (request.journal)
  {
    for (Families::Iterator it(families); it.More(); it.Next())
    {
      if (it.Value().IsEmpty())
        result.history->Remove(it.Key());
      for (const auto& product : it.Value())
        if (!product.IsSame(it.Key()))
          result.history->AddModified(it.Key(), product);
    }
    // Intersections generated from input boundaries are composed through reunion too.
    // BOPAlgo reports these separately from Modified; only final entities can be published.
    ShapeSet final;
    TopExp::MapShapes(result.shape, final);
    for (const auto& operand : partition.operands)
      for (auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
        for (TopExp_Explorer it(operand.original, type); it.More(); it.Next())
          for (const auto& generated : partition.builder->Generated(it.Current()))
            if (final.Contains(generated))
              result.history->AddGenerated(it.Current(), generated);
            else
              for (const auto& successor : transformation->Modified(generated))
                if (final.Contains(successor))
                  result.history->AddGenerated(it.Current(), successor);
  }
  return result;
}

} // namespace

BooleanResult booleanCells(const BooleanRequest& request, std::vector<Rung>* report)
{
  checkCancelled(request.operation, request.cancellable);
  BooleanPhase ladderPhase("ladder");
  const double preparation = booleanTimingEnabled() ? booleanNow() : 0;
  if (request.sources.empty())
    fail(request.operation, invalid_input, "nothing to combine");
  // OCCT's temporary containers freeze the inputs and their parts. Restore their incoming
  // flags after those containers die. The returned bodies have owned roots (see build); a
  // part they share with an input is left as free as the input has it.
  Freedom freedom(request);
  if (request.select == Select::Region)
  {
    int dimension = 0;
    for (const auto* group : {&request.sources, &request.tools})
      for (const auto& shape : *group)
      {
        const int next = dimensionOf(request.operation, shape);
        if (next < 2)
          fail(request.operation,
               unsupported_geometry,
               "regions require solids or sheets, not curves");
        if (dimension && dimension != next)
          fail(request.operation,
               unsupported_geometry,
               "mixed solid/sheet Region imprints are not implemented");
        dimension = next;
      }
  }
  double  tolerance = Precision::Confusion(), diagonal = 0;
  Bnd_Box box;
  for (const auto* group : {&request.sources, &request.tools})
    for (const auto& shape : *group)
    {
      dimensionOf(request.operation, shape);
      tolerance = std::max({tolerance, tubeOf(shape), vertexToleranceOf(shape)});
      BRepBndLib::Add(shape, box, false);
    }
  if (!box.IsVoid() && !box.IsOpen())
    diagonal = std::sqrt(box.SquareExtent());

  // The ladder. Each rung is verified against the operands as the caller gave them.
  ShapeList            lines;
  std::vector<Operand> operands;
  {
    for (const auto& shape : request.sources)
      operands.push_back(operandOf(request, shape, false));
    for (const auto& shape : request.tools)
      operands.push_back(operandOf(request, shape, true));
    BooleanPhase crossings("sweptCrossings");
    lines = sweptCrossings(operands, tolerance::booleanCrossing);
  }
  spherePoles(operands);
  bool   torusPair = operands.size() == 2;
  gp_Dir torusAxes[2];
  size_t torusIndex = 0;
  if (torusPair)
    for (const auto& operand : operands)
    {
      TopExp_Explorer face(operand.original, TopAbs_FACE);
      if (!face.More()
          || BRepAdaptor_Surface(TopoDS::Face(face.Current())).GetType() != GeomAbs_Torus)
      {
        torusPair = false;
        continue;
      }
      torusAxes[torusIndex++] =
        BRepAdaptor_Surface(TopoDS::Face(face.Current())).Torus().Axis().Direction();
      face.Next();
      if (face.More())
        torusPair = false;
    }
  if (torusPair && !torusAxes[0].IsParallel(torusAxes[1], BRepExt::tolerance::angularResolution))
    divideSeams(operands);
  // Check 4 does not yet bound lower-dimensional contacts against faces. Do not let
  // a caller's fuzzy value (including a command retry) erase contacts above the ceiling.
  const double      ceiling     = tolerance::booleanCeiling;
  const double      callerFuzzy = std::min(request.fuzzy, ceiling);
  std::vector<Step> ladder{{"plain", callerFuzzy}};
  const double      fuzzy =
    std::min({kFuzzyFactor * tolerance, tolerance::booleanFuzzyRelative * diagonal, ceiling});
  if (fuzzy > callerFuzzy)
    ladder.push_back({"fuzzy", fuzzy});

  if (booleanTimingEnabled())
    booleanPhase("preparation", since(preparation));
  std::vector<Rung> rungs;
  BooleanResult     result;
  bool              made = false;
  for (const auto& step : ladder)
  {
    checkCancelled(request.operation, request.cancellable);
    rungs.emplace_back();
    auto& rung = rungs.back();
    rung.name  = step.name;
    rung.fuzzy = step.fuzzy;
    BooleanPhase rungPhase("rung");

    // Declare the destruction timer before the partition so it can include its lifetime end.
    struct Destruction
    {
      double start = 0;

      ~Destruction()
      {
        if (start)
          booleanPhase("destruction", since(start));
      }
    } destruction;

    Partition partition;
    auto      start = Clock::now();
    try
    {
      try
      {
        fuse(partition, request, step, lines, &operands);
        rung.fuseMs    = since(start);
        rung.cells     = static_cast<int>(partition.cells.size());
        rung.tolerance = partition.tolerance;
        start          = Clock::now();
        verify(partition, request, rung);
        rung.checkMs  = since(start);
        start         = Clock::now();
        result        = build(partition, request, rung);
        rung.buildMs  = since(start);
        rung.accepted = made = true;
      }
      catch (const Standard_Failure& error)
      {
        refuse(std::string("OpenCascade failed (") + error.what() + ")");
      }
    }
    catch (const Refused& refused)
    {
      rung.refusal = refused.why;
      if (rung.fuseMs == 0)
        rung.fuseMs = since(start);
      else if (rung.checkMs == 0)
        rung.checkMs = since(start);
    }
    operands = std::move(partition.operands);
    if (booleanTimingEnabled())
      destruction.start = booleanNow();
    if (made)
      break;
  }
  if (report)
    *report = rungs;
  if (!made)
  {
    std::string why = rungs.front().refusal;
    for (size_t i = 1; i < rungs.size(); ++i)
      if (rungs[i].refusal != rungs.front().refusal)
        why += "; " + rungs[i].name + ": " + rungs[i].refusal;
    why += ". Move or rotate the tool to change the intersection, or make the overlapping feature "
           "thicker. Changing fuzzy tolerance cannot replace a verified partition.";
    fail(request.operation, no_solution, why);
  }
  result.rungs = rungs;
  return result;
}

} // namespace BRepExt
