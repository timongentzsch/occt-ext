// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
// Draft of walls that OpenCascade's draft does not take (BRepOffsetAPI_DraftAngle inclines
// planes, cylinders and cones only): a wall along a spline, or along any curve.
//
// The drafted wall is the ruled isocline surface over the wall's neutral curve: at the height
// z above the neutral plane it is the curve moved sideways, along its own normal in that
// plane, by -z tan(angle). Every tangent plane of it makes the angle with the pull direction.
// The wall is replaced by that surface the way push_face.cpp moves a face: what lies between
// the two is a solid of its own, taken from the body where the new surface lies inside it
// and joined to it where it lies outside.
//
// What is taken stands a little out of the body, past the wall and past its upper and lower
// ends, and what is joined stands a little into it, so that no face of either lies in the wall
// itself. Their ends along the wall, and the end of a joined one at the wall's upper or lower
// edge, do lie in the neighbouring faces if those are flat and square to the wall there; the
// next wall of the same curve, tangent to this one, shares that end. Where an end does not
// lie in a neighbour it would be left as a step in the result, and the draft is refused.
//
// For walls along the pull direction between two levels: the side faces of an extrusion.

#include <BRepSkinExt_DraftIsocline.hxx>
#include <BRepOffsetExt_Solid.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_OffsetCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <cmath>

namespace BRepExt
{

namespace
{

const char* const op = "draftFaces";

/// The curves of the drafted wall and of the walls beside it are fitted to this: the
/// distance at which OpenCascade takes two points for one.
constexpr double isoclineFit = BRepExt::tolerance::length;
/// The faces of a solid that is taken or joined are sewn at ten times that fit.
constexpr double isoclineSew = BRepExt::tolerance::sewing;
/// How far what is taken stands out of the body and what is joined into it, in units of the
/// caller's tolerance: far above the Boolean's fuzz (1e-7), small against any wall.
constexpr double isoclineOverlap = 1e3;
/// How far (radians) a wall's normal may lean from square to the pull direction and the wall
/// still count as running along it, and how far the angle that is measured on the result may
/// be from the angle asked: the normals of a surface ruled between curves fitted to 1e-7.
constexpr double isoclineAngle = 1e-5;

/// One wall as the draft needs it: a curve at the wall's lower level, the levels, and the
/// side of the curve that the wall faces.
struct Wall
{
  TopoDS_Face        face;
  Handle(Geom_Curve) curve;             // at `low`, in a plane across the pull direction
  double             low = 0, high = 0; // heights over the neutral plane
  double             side   = 1;        // the wall faces along side * (tangent x pull)
  bool               closed = false;
};

double heightOf(const gp_Pnt& point, const gp_Pln& neutral, const gp_Dir& pull)
{
  return gp_Vec(neutral.Location(), point).Dot(gp_Vec(pull));
}

/// The wall, or a failure that names the face: it has to run along the pull direction, from
/// one level to another, between straight ends.
Wall wallOf(const TopoDS_Face& face, const gp_Pln& neutral, const gp_Dir& pull, double tolerance)
{
  Wall wall;
  wall.face = face;
  const BRepAdaptor_Surface surface(face);
  gp_Dir                    outward;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
    {
      const double u = surface.FirstUParameter()
                       + (surface.LastUParameter() - surface.FirstUParameter()) * (i + 0.5) / 5;
      const double v = surface.FirstVParameter()
                       + (surface.LastVParameter() - surface.FirstVParameter()) * (j + 0.5) / 5;
      gp_Pnt at;
      gp_Dir normal;
      if (!outwardNormal(face, u, v, at, normal) || std::abs(normal.Dot(pull)) > isoclineAngle)
        fail(op,
             unsupported_geometry,
             "a curved face can be drafted only where it runs along the pull direction (the wall "
             "of an extrusion)",
             {face});
    }
  int wires = 0;
  for (TopExp_Explorer all(face, TopAbs_WIRE); all.More(); all.Next())
    wires++;
  // Its edges: along the pull direction, or at one of two levels.
  std::vector<TopoDS_Edge> across;
  std::vector<double>      levels;
  bool                     first = true;
  for (TopExp_Explorer all(face, TopAbs_EDGE); all.More(); all.Next())
  {
    const TopoDS_Edge edge = TopoDS::Edge(all.Current());
    if (BRep_Tool::Degenerated(edge))
      continue;
    const BRepAdaptor_Curve curve(edge);
    double                  least = 0, most = 0;
    for (int i = 0; i <= 8; i++)
    {
      const double z =
        heightOf(curve.Value(curve.FirstParameter()
                             + (curve.LastParameter() - curve.FirstParameter()) * i / 8),
                 neutral,
                 pull);
      least = i == 0 ? z : std::min(least, z);
      most  = i == 0 ? z : std::max(most, z);
    }
    if (most - least <= tolerance)
    {
      across.push_back(edge);
      levels.push_back(0.5 * (least + most));
    }
    else if (curve.GetType() != GeomAbs_Line
             || std::abs(curve.Line().Direction().Dot(pull)) < 1.0 - isoclineAngle)
    {
      wires = 0;
    }
    wall.low  = first ? least : std::min(wall.low, least);
    wall.high = first ? most : std::max(wall.high, most);
    first     = false;
  }
  for (double level : levels)
    if (std::abs(level - wall.low) > tolerance && std::abs(level - wall.high) > tolerance)
      wires = 0;
  if (wires != 1 || wall.high - wall.low <= tolerance)
    fail(op,
         unsupported_geometry,
         "a curved face can be drafted only if it is a wall between two levels across the pull "
         "direction, without holes",
         {face});

  // The curve at the lower level, of one piece.
  GeomConvert_CompCurveToBSplineCurve joined;
  std::vector<TopoDS_Edge>            lower;
  for (size_t i = 0; i < across.size(); i++)
    if (std::abs(levels[i] - wall.low) <= tolerance)
      lower.push_back(across[i]);
  // In the order they follow each other: an edge is added where it fits an end.
  std::vector<bool> used(lower.size(), false);
  for (size_t added = 0; added < lower.size(); added++)
  {
    bool fitted = false;
    for (size_t i = 0; i < lower.size() && !fitted; i++)
    {
      if (used[i])
        continue;
      double                   from, to;
      const Handle(Geom_Curve) piece = BRep_Tool::Curve(lower[i], from, to);
      if (piece.IsNull())
        continue;
      const Handle(Geom_BSplineCurve) spline =
        GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(piece, from, to));
      fitted =
        joined.Add(spline, 10 * tolerance, true) || joined.Add(spline, 10 * tolerance, false);
      used[i] = fitted;
    }
    if (!fitted)
      fail(op, unsupported_geometry, "the lower edge of a wall is not one curve", {face});
  }
  if (lower.empty())
    fail(op, unsupported_geometry, "a wall has no edge at its lower level", {face});
  wall.curve  = joined.BSplineCurve();
  wall.closed = wall.curve->Value(wall.curve->FirstParameter())
                  .Distance(wall.curve->Value(wall.curve->LastParameter()))
                <= 10 * tolerance;

  // Which side of the curve the wall faces: where the curve's own normal agrees with the face's.
  const double mid = 0.5 * (wall.curve->FirstParameter() + wall.curve->LastParameter());
  gp_Pnt       on;
  gp_Vec       tangent;
  wall.curve->D1(mid, on, tangent);
  const gp_Vec across2 = tangent.Crossed(gp_Vec(pull));
  // The face's normal there: at the point of the face nearest to the curve's middle, lifted to mid
  // height.
  double best = -1;
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
    {
      const double u = surface.FirstUParameter()
                       + (surface.LastUParameter() - surface.FirstUParameter()) * (i + 0.5) / 9;
      const double v = surface.FirstVParameter()
                       + (surface.LastVParameter() - surface.FirstVParameter()) * (j + 0.5) / 9;
      gp_Pnt at;
      gp_Dir normal;
      if (!outwardNormal(face, u, v, at, normal))
        continue;
      const gp_Vec offset(on, at);
      const double sideways = (offset - gp_Vec(pull) * offset.Dot(gp_Vec(pull))).Magnitude();
      if (best < 0 || sideways < best)
      {
        best    = sideways;
        outward = normal;
      }
    }
  wall.side = gp_Vec(outward).Dot(across2) >= 0 ? 1.0 : -1.0;
  return wall;
}

/// The wall's curve at the height `z`, moved sideways out of the wall by `shift`: fitted.
Handle(Geom_BSplineCurve) shifted(const Wall& wall, const gp_Dir& pull, double z, double shift)
{
  Handle(Geom_Curve) at =
    Handle(Geom_Curve)::DownCast(wall.curve->Translated(gp_Vec(pull) * (z - wall.low)));
  if (std::abs(shift) <= gp::Resolution())
    return GeomConvert::CurveToBSplineCurve(at);
  // Geom_OffsetCurve moves along tangent x direction.
  const Handle(Geom_OffsetCurve) moved = new Geom_OffsetCurve(at, wall.side * shift, pull);
  for (const GeomAbs_Shape order : {GeomAbs_C2, GeomAbs_C1, GeomAbs_C0})
  {
    GeomConvert_ApproxCurve fit(moved, isoclineFit, order, 400, 9);
    if (fit.HasResult() && fit.MaxError() <= isoclineFit)
      return fit.Curve();
  }
  return {};
}

/// The least of 1 - shift * curvature along the wall's curve: at zero and below the curve
/// moved by `shift` folds over itself.
double roomFor(const Wall& wall, const gp_Dir& pull, double shift)
{
  double       least = 1.0;
  const double from = wall.curve->FirstParameter(), to = wall.curve->LastParameter();
  for (int i = 0; i <= 256; i++)
  {
    gp_Pnt at;
    gp_Vec d1, d2;
    wall.curve->D2(from + (to - from) * i / 256, at, d1, d2);
    const double speed = d1.Magnitude();
    if (speed <= gp::Resolution())
      continue;
    const gp_Vec normal = d1.Crossed(gp_Vec(pull)) / speed * wall.side;
    // The curvature towards the side the wall faces.
    const double towards = d2.Dot(normal) / (speed * speed);
    least                = std::min(least, 1.0 - shift * towards);
  }
  return least;
}

struct Ruled
{
  Handle(Geom_Surface) surface;
  bool                 alongU = true; // the curve runs along u, the height along v
  double               from = 0, to = 0;

  /// The (u, v) limits of the piece between the heights a and b of a surface that spans z0 to z1.
  TopoDS_Face piece(double z0, double z1, double a, double b) const
  {
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    const double low = alongU ? v0 : u0, high = alongU ? v1 : u1;
    const double start = low + (high - low) * (a - z0) / (z1 - z0),
                 end   = low + (high - low) * (b - z0) / (z1 - z0);
    return alongU ? BRepBuilderAPI_MakeFace(surface, u0, u1, start, end, isoclineFit).Face()
                  : BRepBuilderAPI_MakeFace(surface, start, end, v0, v1, isoclineFit).Face();
  }

  Handle(Geom_Curve) level(double z0, double z1, double z) const
  {
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    const double low = alongU ? v0 : u0, high = alongU ? v1 : u1;
    const double at = low + (high - low) * (z - z0) / (z1 - z0);
    return alongU ? surface->VIso(at) : surface->UIso(at);
  }
};

/// The surface ruled between two curves, and which way its parameters run.
Ruled ruledBetween(const Handle(Geom_BSplineCurve)& lower, const Handle(Geom_BSplineCurve)& upper)
{
  Ruled                  ruled;
  GeomFill_BSplineCurves fill(lower, upper, GeomFill_StretchStyle);
  ruled.surface = fill.Surface();
  if (ruled.surface.IsNull())
    return ruled;
  double u0, u1, v0, v1;
  ruled.surface->Bounds(u0, u1, v0, v1);
  // The lower curve is an isoparametric line of it: find out which.
  const gp_Pnt far    = lower->Value(0.5 * (lower->FirstParameter() + lower->LastParameter()));
  const gp_Pnt alongU = ruled.surface->Value(0.5 * (u0 + u1), v0),
               alongV = ruled.surface->Value(u0, 0.5 * (v0 + v1));
  ruled.alongU        = alongU.Distance(far) <= alongV.Distance(far);
  return ruled;
}

/// A solid between two walls ruled over the same curve, from the height a to b.
TopoDS_Shape slab(const Ruled& inner,
                  const Ruled& outer,
                  double       z0,
                  double       z1,
                  double       a,
                  double       b,
                  bool         closed)
{
  BRepBuilderAPI_Sewing sewing(isoclineSew);
  sewing.Add(inner.piece(z0, z1, a, b));
  sewing.Add(outer.piece(z0, z1, a, b));
  std::vector<gp_Pnt> corners[2];
  for (double z : {a, b})
  {
    const Handle(Geom_Curve) in = inner.level(z0, z1, z), out = outer.level(z0, z1, z);
    const TopoDS_Edge        first = BRepBuilderAPI_MakeEdge(in).Edge(),
                      second       = BRepBuilderAPI_MakeEdge(out).Edge();
    sewing.Add(BRepFill::Face(first, second));
    for (int end = 0; end < 2; end++)
    {
      corners[end].push_back(in->Value(end == 0 ? in->FirstParameter() : in->LastParameter()));
      corners[end].push_back(out->Value(end == 0 ? out->FirstParameter() : out->LastParameter()));
    }
  }
  if (!closed)
    for (int end = 0; end < 2; end++)
    {
      // In at a, out at a, out at b, in at b: round the end.
      BRepBuilderAPI_MakePolygon polygon(corners[end][0],
                                         corners[end][1],
                                         corners[end][3],
                                         corners[end][2],
                                         true);
      if (!polygon.IsDone())
        return {};
      BRepBuilderAPI_MakeFace face(polygon.Wire(), true);
      if (!face.IsDone())
        return {};
      sewing.Add(face.Face());
    }
  sewing.Perform();
  const TopoDS_Shape sewn = sewing.SewedShape();
  if (sewn.IsNull() || sewn.ShapeType() != TopAbs_SHELL)
    return {};
  TopoDS_Solid solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(sewn)).Solid();
  if (!BRepLib::OrientClosedSolid(solid) || !BRepCheck_Analyzer(solid).IsValid()
      || !(estimatedVolume(solid) > 0.0))
    return {};
  return solid;
}

} // namespace

TopoDS_Shape draftWalls(const TopoDS_Shape&              solid,
                        const std::vector<TopoDS_Shape>& faces,
                        double                           angle,
                        const gp_Dir&                    pull,
                        const gp_Pln&                    neutral,
                        double                           tolerance,
                        Lineage&                         lineage)
{
  if (std::abs(neutral.Axis().Direction().Dot(pull)) < 1.0 - isoclineAngle)
    fail(op,
         unsupported_geometry,
         "a curved face is drafted about a neutral plane across the pull direction only",
         faces);
  const double overlap = isoclineOverlap * tolerance;
  const double slope   = std::tan(angle);

  TopoDS_Compound tools;
  BRep_Builder().MakeCompound(tools);
  NCollection_List<TopoDS_Shape>    taken, joined;
  std::vector<Handle(Geom_Surface)> drafted;
  for (const TopoDS_Shape& shape : faces)
  {
    const Wall wall = wallOf(TopoDS::Face(shape), neutral, pull, tolerance);
    // The wall over the whole height, and the walls just inside and outside it.
    const double z0 = wall.low - overlap, z1 = wall.high + overlap;
    const auto   shiftAt = [&](double z) { return -z * slope; };
    // The walls beside it, far enough that the drafted wall stays between them where it
    // runs on past the neutral plane.
    const double beside = overlap * (1.0 + std::abs(slope));
    // Where the curve is bent towards the side it is moved to more than it is moved, it folds.
    if (std::min(roomFor(wall, pull, shiftAt(z0)), roomFor(wall, pull, shiftAt(z1)))
        <= collapseSlack)
      fail(op,
           no_solution,
           "the taper moves a wall further than its radius of curvature; reduce the angle",
           {shape});
    const Handle(Geom_BSplineCurve) lowR    = shifted(wall, pull, z0, shiftAt(z0)),
                                    highR   = shifted(wall, pull, z1, shiftAt(z1));
    const Handle(Geom_BSplineCurve) lowIn   = shifted(wall, pull, z0, -beside),
                                    highIn  = shifted(wall, pull, z1, -beside);
    const Handle(Geom_BSplineCurve) lowOut  = shifted(wall, pull, z0, beside),
                                    highOut = shifted(wall, pull, z1, beside);
    if (lowR.IsNull() || highR.IsNull() || lowIn.IsNull() || highIn.IsNull() || lowOut.IsNull()
        || highOut.IsNull())
      fail(op, no_solution, "the curve of a wall could not be moved sideways to 1e-7", {shape});
    const Ruled isocline = ruledBetween(lowR, highR), inside = ruledBetween(lowIn, highIn),
                outside = ruledBetween(lowOut, highOut);
    if (isocline.surface.IsNull() || inside.surface.IsNull() || outside.surface.IsNull())
      fail(op, no_solution, "the drafted wall could not be built", {shape});
    drafted.push_back(isocline.surface);

    // Below and above the neutral plane the wall leans to opposite sides.
    struct Part
    {
      double from, to;
    };

    std::vector<Part> parts;
    if (wall.low < -tolerance)
      parts.push_back({wall.low, std::min(wall.high, 0.0)});
    if (wall.high > tolerance)
      parts.push_back({std::max(wall.low, 0.0), wall.high});
    for (const Part& part : parts)
    {
      const double middle   = 0.5 * (part.from + part.to);
      const bool   outwards = shiftAt(middle) > 0;
      // What is joined ends at the wall's own end; what is taken, and either at the
      // neutral plane, runs on by the overlap.
      const bool endsLow      = std::abs(part.from - wall.low) <= tolerance,
                 endsHigh     = std::abs(part.to - wall.high) <= tolerance;
      const double       a    = outwards && endsLow ? part.from : part.from - overlap;
      const double       b    = outwards && endsHigh ? part.to : part.to + overlap;
      const TopoDS_Shape tool = outwards ? slab(inside, isocline, z0, z1, a, b, wall.closed)
                                         : slab(isocline, outside, z0, z1, a, b, wall.closed);
      if (tool.IsNull())
        fail(op, no_solution, "the solid between a wall and its draft could not be built", {shape});
      (outwards ? joined : taken).Append(tool);
      BRep_Builder().Add(tools, tool);
    }
  }

  // Take, then join.
  TopoDS_Shape                   result = solid;
  Lineage                        swept(tools);
  NCollection_List<TopoDS_Shape> body;
  body.Append(result);
  if (!taken.IsEmpty())
  {
    BRepAlgoAPI_Cut cut;
    cut.SetArguments(body);
    cut.SetTools(taken);
    cut.SetNonDestructive(true);
    cut.Build();
    if (!cut.IsDone() || cut.HasErrors() || singleSolid(cut.Shape()).IsNull())
      fail(op, no_solution, "what the taper removes could not be taken from the body", faces);
    lineage.through(cut);
    swept.through(cut);
    result = cut.Shape();
  }
  if (!joined.IsEmpty())
  {
    body.Clear();
    body.Append(result);
    BRepAlgoAPI_Fuse fuse;
    fuse.SetArguments(body);
    fuse.SetTools(joined);
    fuse.SetNonDestructive(true);
    fuse.Build();
    if (!fuse.IsDone() || fuse.HasErrors() || singleSolid(fuse.Shape()).IsNull())
      fail(op, no_solution, "what the taper adds could not be joined to the body", faces);
    lineage.through(fuse);
    swept.through(fuse);
    result = fuse.Shape();
  }
  // One face again where what was joined lies in a neighbour; the edges that the body had are kept.
  ShapeUpgrade_UnifySameDomain                           unify(result, true, true, false);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> keep;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> ofWalls;
  for (const TopoDS_Shape& face : faces)
    for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
      ofWalls.Add(edges.Current());
  for (TopExp_Explorer edges(solid, TopAbs_EDGE); edges.More(); edges.Next())
    if (!ofWalls.Contains(edges.Current()))
      keep.Add(edges.Current());
  unify.KeepShapes(keep);
  unify.Build();
  result = singleSolid(unify.Shape());
  if (result.IsNull() || !BRepCheck_Analyzer(result).IsValid())
    fail(op, no_solution, "the drafted body is not one valid solid", faces);
  lineage.throughHistory(unify.History());
  swept.throughHistory(unify.History());
  lineage.within(result);
  swept.within(result);

  // Each face of the result that came with a solid of ours is the drafted wall, or lies in a
  // face that the body had; anything else is a step that the draft would leave behind.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> known;
  for (const TopoDS_Shape& face : lineage.all())
    known.Add(face);
  std::vector<std::vector<TopoDS_Shape>> walls(faces.size());
  for (const TopoDS_Shape& face : swept.all())
  {
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(TopoDS::Face(face));
    size_t                     which   = drafted.size();
    for (size_t i = 0; i < drafted.size(); i++)
      if (drafted[i] == surface)
        which = i;
    if (which < drafted.size())
    {
      bool seen = false;
      for (const TopoDS_Shape& have : walls[which])
        seen = seen || have.IsSame(face);
      if (!seen)
        walls[which].push_back(face);
    }
    else if (!known.Contains(face))
    {
      fail(op,
           unsupported_geometry,
           "the draft would leave a step: the faces beside a drafted wall must be flat and square "
           "to it where they meet, or be drafted with it and tangent to it",
           faces);
    }
  }
  // Apart from the Boolean's word: the angle of the drafted walls to the pull direction.
  for (size_t i = 0; i < faces.size(); i++)
  {
    if (walls[i].empty())
      fail(op, no_solution, "a drafted wall is gone from the body", {faces[i]});
    for (const TopoDS_Shape& face : walls[i])
    {
      const BRepAdaptor_Surface surface(TopoDS::Face(face));
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
        {
          gp_Pnt       point;
          gp_Dir       normal;
          const double u = surface.FirstUParameter()
                           + (surface.LastUParameter() - surface.FirstUParameter()) * (a + 0.5) / 3;
          const double v = surface.FirstVParameter()
                           + (surface.LastVParameter() - surface.FirstVParameter()) * (b + 0.5) / 3;
          if (outwardNormal(TopoDS::Face(face), u, v, point, normal)
              && std::abs(std::abs(std::asin(std::clamp(normal.Dot(pull), -1.0, 1.0)))
                          - std::abs(angle))
                   > isoclineAngle)
            fail(op, no_solution, "a drafted wall does not have the angle", {faces[i]});
        }
    }
    lineage.set(faces[i], walls[i]);
  }
  return result;
}

} // namespace BRepExt
