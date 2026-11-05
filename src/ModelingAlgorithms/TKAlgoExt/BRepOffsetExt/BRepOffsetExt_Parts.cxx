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

#include <BRepOffsetExt_Parts.hxx>

#include <BRepSurfaceExt_FaceSeparation.hxx>

#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <GeomLib.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Adaptor3d_CurveOnSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom2d_Curve.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepBndLib.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <deque>
#include <BRepLProp_SLProps.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS_Compound.hxx>
#include <ChFi3d.hxx>
#include <IntTools_Context.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

#include <cmath>
#include <limits>

namespace BRepExt
{

namespace
{

//! The (u, v) inside `face` on a grid of `perSide` squared, cell centres; along lines across it
//! where the grid misses it (a ring, a slot).
std::vector<gp_Pnt2d> parametersInside(const TopoDS_Face& face, int perSide)
{
  const BRepAdaptor_Surface surface(face);
  BRepTopAdaptor_FClass2d   inside(face, Precision::Confusion());
  const double              u0 = surface.FirstUParameter(), u1 = surface.LastUParameter();
  const double              v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
  std::vector<gp_Pnt2d>     out;
  for (int i = 0; i < perSide; i++)
    for (int j = 0; j < perSide; j++)
    {
      const gp_Pnt2d uv(u0 + (u1 - u0) * (i + 0.5) / perSide, v0 + (v1 - v0) * (j + 0.5) / perSide);
      if (inside.Perform(uv) == TopAbs_IN)
        out.push_back(uv);
    }
  if (out.empty())
    out = parametersAcross(face, perSide, new IntTools_Context);
  return out;
}

NCollection_List<TopoDS_Shape> listOf(const std::vector<TopoDS_Shape>& faces)
{
  NCollection_List<TopoDS_Shape> list;
  for (const TopoDS_Shape& face : faces)
    list.Append(face);
  return list;
}

//! A distance query must use the trimmed surface, not an inherited spatial
//! edge curve which may lie anywhere inside its input tolerance tube.
TopoDS_Shape surfaceTrims(const TopoDS_Shape& input)
{
  TopoDS_Compound out;
  BRep_Builder    builder;
  builder.MakeCompound(out);
  for (TopExp_Explorer ff(input, TopAbs_FACE); ff.More(); ff.Next())
  {
    BRepBuilderAPI_Copy copy(ff.Current(), true, false);
    const auto          face    = TopoDS::Face(copy.Shape());
    const auto          support = BRep_Tool::Surface(face);
    for (TopExp_Explorer ee(face, TopAbs_EDGE); ee.More(); ee.Next())
    {
      const auto edge = TopoDS::Edge(ee.Current());
      if (BRep_Tool::Degenerated(edge))
        continue;
      double     first, last;
      const auto on = BRep_Tool::CurveOnSurface(edge, face, first, last);
      if (on.IsNull())
        throw Standard_ConstructionError("missing surface trim");
      TopoDS_Vertex a, b;
      TopExp::Vertices(edge, a, b);
      const auto move = [&](const TopoDS_Vertex& vertex, double t) {
        if (vertex.IsNull())
          return;
        const auto uv = on->Value(t);
        builder.UpdateVertex(vertex, support->Value(uv.X(), uv.Y()), Precision::Confusion());
      };
      move(a, first);
      move(b, last);
      const Handle(Geom2dAdaptor_Curve) curve2d = new Geom2dAdaptor_Curve(on, first, last);
      const Handle(GeomAdaptor_Surface) surface = new GeomAdaptor_Surface(support);
      Adaptor3d_CurveOnSurface          trim(curve2d, surface);
      Handle(Geom_Curve)                curve3d;
      double                            maximum = 0, average = 0;
      GeomLib::BuildCurve3d(0.1 * Precision::Confusion(),
                            trim,
                            first,
                            last,
                            curve3d,
                            maximum,
                            average,
                            GeomAbs_C1,
                            14,
                            1000);
      if (curve3d.IsNull() || !std::isfinite(maximum) || maximum > 0.1 * Precision::Confusion())
        throw Standard_ConstructionError("surface trim approximation unresolved");
      builder.UpdateEdge(edge, curve3d, Precision::Confusion());
    }
    builder.Add(out, face);
  }
  return out;
}

} // namespace

OffsetFate classifyForOffset(const TopoDS_Shape& body, double distance)
{
  OffsetFate fate;
  for (TopExp_Explorer faces(body, TopAbs_FACE); faces.More(); faces.Next())
  {
    const TopoDS_Face&        face = TopoDS::Face(faces.Current());
    const BRepAdaptor_Surface surface(face);
    if (surface.GetType() == GeomAbs_Plane)
      continue;
    // Along the surface's own normal the offset is this far.
    const double along = face.Orientation() == TopAbs_REVERSED ? -distance : distance;
    // Elementary radii are exact data, not fitted blend estimates. Compare
    // physical focal distances directly to avoid the relative collapse allowance.
    if (surface.GetType() == GeomAbs_Sphere || surface.GetType() == GeomAbs_Cylinder)
    {
      const double radius = surface.GetType() == GeomAbs_Sphere ? surface.Sphere().Radius()
                                                                : surface.Cylinder().Radius();
      if (-along >= radius)
        fate.vanish.push_back(face);
      continue;
    }
    if (surface.GetType() == GeomAbs_Cone)
    {
      const auto   cone = surface.Cone();
      const double r0   = cone.RefRadius() + surface.FirstVParameter() * std::sin(cone.SemiAngle());
      const double r1   = cone.RefRadius() + surface.LastVParameter() * std::sin(cone.SemiAngle());
      const double minimum = r0 * r1 <= 0 ? 0 : std::min(std::abs(r0), std::abs(r1));
      const double maximum = std::max(std::abs(r0), std::abs(r1));
      const double reached = -along * std::cos(cone.SemiAngle());
      if (reached >= maximum)
        fate.vanish.push_back(face);
      else if (reached >= minimum)
        fate.partial.push_back(face);
      continue;
    }
    if (surface.GetType() == GeomAbs_Torus
        && surface.Torus().MajorRadius() > surface.Torus().MinorRadius())
    {
      const double major = surface.Torus().MajorRadius(), minor = surface.Torus().MinorRadius();
      if (-along >= minor)
        fate.vanish.push_back(face);
      else if (along > 0)
      {
        const double v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
        double       low  = std::min(std::cos(v0), std::cos(v1));
        double       high = std::max(std::cos(v0), std::cos(v1));
        const double turn = 2 * M_PI;
        if (v1 - v0 >= turn)
        {
          low  = -1;
          high = 1;
        }
        else
        {
          if (std::ceil(v0 / turn) <= std::floor(v1 / turn))
            high = 1;
          if (std::ceil((v0 - M_PI) / turn) <= std::floor((v1 - M_PI) / turn))
            low = -1;
        }
        const double threshold = -major / (minor + along);
        if (high <= threshold)
          fate.vanish.push_back(face);
        else if (low <= threshold)
          fate.partial.push_back(face);
      }
      continue;
    }
    const bool elementary = surface.GetType() == GeomAbs_Torus;
    int        collapsed = 0, measured = 0;
    // 9 x 9: a blend is the same all along, and a spline face that is tight in a part
    // smaller than a ninth of its width is left to the checks on the result.
    for (const gp_Pnt2d& uv : parametersInside(face, 9))
    {
      BRepLProp_SLProps at(surface, uv.X(), uv.Y(), 2, Precision::Confusion());
      if (!at.IsCurvatureDefined())
        continue;
      measured++;
      // A principal curvature k becomes k / (1 - along * k): gone where along * k is 1.
      if (std::max(along * at.MaxCurvature(), along * at.MinCurvature())
          >= (elementary ? 1.0 : 1.0 - collapseSlack))
        collapsed++;
    }
    if (collapsed == 0)
      continue;
    (collapsed == measured ? fate.vanish : fate.partial).push_back(face);
  }
  return fate;
}

bool intersectsItself(const TopoDS_Shape& shape)
{
  // Only what the separation pre-pass cannot prove clear goes to the check
  // (BRepSurfaceExt_FaceSeparation.hxx).
  for (const TopoDS_Shape& part : unprovenParts(shape))
  {
    BOPAlgo_ArgumentAnalyzer analyzer;
    analyzer.SetShape1(part);
    analyzer.ArgumentTypeMode() = false;
    analyzer.SmallEdgeMode()    = false;
    analyzer.SelfInterMode()    = true;
    analyzer.Perform();
    if (analyzer.HasFaulty())
      return true;
  }
  return false;
}

std::vector<gp_Pnt> pointsInside(const TopoDS_Face& face, int perSide)
{
  const BRepAdaptor_Surface surface(face);
  std::vector<gp_Pnt>       out;
  for (const gp_Pnt2d& uv : parametersInside(face, perSide))
    out.push_back(surface.Value(uv.X(), uv.Y()));
  return out;
}

Clearance clearance(const std::vector<TopoDS_Shape>& faces,
                    const TopoDS_Shape&              body,
                    const std::vector<TopoDS_Shape>& opened,
                    double                           wall,
                    double                           required)
{
  TopoDS_Compound skin, openings;
  BRep_Builder    builder;
  builder.MakeCompound(skin);
  builder.MakeCompound(openings);
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next())
  {
    bool open = false;
    for (const auto& face : opened)
      open = open || face.IsSame(all.Current());
    builder.Add(open ? openings : skin, all.Current());
  }
  Clearance result;
  result.expected = faces.size();
  if (faces.empty() || !TopExp_Explorer(skin, TopAbs_FACE).More())
  {
    result.resolved = false;
    result.issue    = "empty skin or request";
    return result;
  }
  try
  {
    skin     = TopoDS::Compound(surfaceTrims(skin));
    openings = TopoDS::Compound(surfaceTrims(openings));
  }
  catch (const Standard_Failure& error)
  {
    result.resolved = false;
    result.issue    = std::string("surface trims: ") + error.what();
    return result;
  }
  DistanceTo toSkin(skin), toOpening(openings);

  struct Cell
  {
    double u0, u1, v0, v1;
    int    depth = 0;
  };

  for (const auto& shape : faces)
  {
    TopoDS_Face face;
    try
    {
      const auto      copy = surfaceTrims(shape);
      TopExp_Explorer ff(copy, TopAbs_FACE);
      if (!ff.More())
      {
        result.resolved = false;
        result.issue    = "missing face";
        return result;
      }
      face = TopoDS::Face(ff.Current());
    }
    catch (const Standard_Failure& error)
    {
      result.resolved = false;
      result.issue    = std::string("surface trim: ") + error.what();
      return result;
    }
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    std::deque<Cell> pending{{u0, u1, v0, v1}};
    int              visited = 0;
    const char*      stage   = "initialization";
    try
    {
      while (!pending.empty())
      {
        checkCancelled("clearance");
        if (++visited > 4096)
        {
          result.resolved = false;
          result.issue    = "subdivision work limit";
          return result;
        }
        const auto cell = pending.front();
        pending.pop_front();
        TopoDS_Shape region = face;
        if (cell.depth)
        {
          stage                = "trim rectangle";
          const auto rectangle = BRepBuilderAPI_MakeFace(BRep_Tool::Surface(face),
                                                         cell.u0,
                                                         cell.u1,
                                                         cell.v0,
                                                         cell.v1,
                                                         Precision::Confusion())
                                   .Face();
          stage                = "trim intersection";
          BRepAlgoAPI_Common common;
          common.SetArguments(listOf({face}));
          common.SetTools(listOf({rectangle}));
          common.SetNonDestructive(true);
          common.Build();
          if (!common.IsDone() || common.HasErrors())
          {
            result.resolved = false;
            result.issue    = "trim intersection failed";
            return result;
          }
          region = common.Shape();
          if (!TopExp_Explorer(region, TopAbs_FACE).More())
            continue;
        }
        stage = "skin distance";
        std::optional<double> distance;
        try
        {
          distance = toSkin.from(region);
        }
        catch (const Standard_Failure&)
        {
        }
        if (distance)
          *distance = std::max(0.0, *distance - 0.2 * Precision::Confusion());
        if (distance && std::isfinite(*distance) && (opened.empty() || *distance >= required))
        {
          result.minimum = std::min(result.minimum, *distance);
          continue;
        }
        // Distance to a closed set is 1-Lipschitz. A box containing the entire
        // trimmed region bounds its distance above by d(center, openings)+radius.
        stage = "region box";
        Bnd_Box box;
        BRepBndLib::AddOptimal(region, box, false, true);
        if (box.IsVoid() || box.IsOpen())
        {
          result.resolved = false;
          result.issue    = "unbounded region";
          return result;
        }
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        const gp_Pnt center((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
        const double radius = std::nextafter(std::hypot(std::hypot(x1 - x0, y1 - y0), z1 - z0) / 2,
                                             std::numeric_limits<double>::infinity());
        // A point distance supplies a conservative lower bound even when the
        // whole-region extremum solver cannot handle this trimmed support.
        if (!distance || !std::isfinite(*distance))
        {
          stage               = "skin point distance";
          const auto atCenter = toSkin.from(center);
          if (atCenter && std::isfinite(*atCenter))
          {
            const double lower = std::nextafter(*atCenter - radius - Precision::Confusion(),
                                                -std::numeric_limits<double>::infinity());
            if (lower >= required)
            {
              result.minimum = std::min(result.minimum, lower);
              continue;
            }
          }
        }
        if (!opened.empty())
        {
          stage             = "opening distance";
          const auto beside = toOpening.from(center);
          if (beside && std::isfinite(*beside)
              && *beside + radius + 1.1 * Precision::Confusion() < wall)
            continue;
        }
        const double um = cell.u0 + (cell.u1 - cell.u0) / 2, vm = cell.v0 + (cell.v1 - cell.v0) / 2;
        if (cell.depth >= 32 || !(um > cell.u0 && um < cell.u1 && vm > cell.v0 && vm < cell.v1))
        {
          result.resolved = false;
          result.issue    = "subdivision precision limit";
          return result;
        }
        Cell lo = cell, hi = cell;
        ++lo.depth;
        ++hi.depth;
        if ((cell.u1 - cell.u0) / (u1 - u0) >= (cell.v1 - cell.v0) / (v1 - v0))
        {
          lo.u1 = um;
          hi.u0 = um;
        }
        else
        {
          lo.v1 = vm;
          hi.v0 = vm;
        }
        pending.push_back(lo);
        pending.push_back(hi);
      }
    }
    catch (const Standard_Failure& error)
    {
      result.resolved = false;
      result.issue    = std::string(stage) + ": " + error.ExceptionType() + ": " + error.what();
      return result;
    }
    ++result.covered;
  }
  return result;
}

bool tangentAlong(const TopoDS_Edge& edge, const TopoDS_Face& a, const TopoDS_Face& b)
{
  return !BRep_Tool::Degenerated(edge) && ChFi3d::IsTangentFaces(edge, a, b);
}

Lineage::Lineage(const TopoDS_Shape& input)
{
  TopExp::MapShapes(input, TopAbs_FACE, input_);
  for (int i = 1; i <= input_.Extent(); i++)
    now_.push_back({input_.FindKey(i)});
}

void Lineage::throughHistory(const Handle(BRepTools_History)& history)
{
  step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
    const NCollection_List<TopoDS_Shape>& next = history->Modified(face);
    if (next.IsEmpty() && !history->IsRemoved(face))
      out.push_back(face);
    for (const TopoDS_Shape& successor : next)
      if (successor.ShapeType() == TopAbs_FACE)
        out.push_back(successor);
  });
}

void Lineage::set(const TopoDS_Shape& face, std::vector<TopoDS_Shape> now)
{
  const int index = input_.FindIndex(face);
  if (index != 0)
    now_[index - 1] = std::move(now);
}

void Lineage::then(const Lineage& next)
{
  step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
    if (!next.input_.Contains(face))
    {
      out.push_back(face);
      return;
    }
    for (const TopoDS_Shape& successor : next.of(face))
      out.push_back(successor);
  });
}

void Lineage::within(const TopoDS_Shape& result)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> faces;
  TopExp::MapShapes(result, TopAbs_FACE, faces);
  step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
    const int index = faces.FindIndex(face);
    if (index == 0)
      return;
    for (const TopoDS_Shape& seen : out)
      if (seen.IsSame(face))
        return;
    out.push_back(faces.FindKey(index));
  });
}

const std::vector<TopoDS_Shape>& Lineage::of(const TopoDS_Shape& face) const
{
  static const std::vector<TopoDS_Shape> none;
  const int                              index = input_.FindIndex(face);
  return index == 0 ? none : now_[index - 1];
}

std::vector<TopoDS_Shape> Lineage::all() const
{
  std::vector<TopoDS_Shape> out;
  for (const auto& faces : now_)
    out.insert(out.end(), faces.begin(), faces.end());
  return out;
}

void Lineage::modified(FaceHistory& history) const
{
  for (int i = 1; i <= input_.Extent(); i++)
  {
    const TopoDS_Shape& face = input_.FindKey(i);
    const auto&         now  = now_[i - 1];
    if (!history.tracks(face) || (now.size() == 1 && now[0].IsSame(face)))
      continue;
    if (now.empty())
      history.deleted(face);
    else
      history.modified(face, listOf(now));
  }
}

void Lineage::generated(FaceHistory& history) const
{
  for (int i = 1; i <= input_.Extent(); i++)
    if (history.tracks(input_.FindKey(i)))
      history.generated(input_.FindKey(i), listOf(now_[i - 1]));
}

} // namespace BRepExt
