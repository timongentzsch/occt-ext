// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepFilletExt_Blend.hxx>
#include <Standard_Failure.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRep_Builder.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <NCollection_Array1.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <ElSLib.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomProjLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <Precision.hxx>
#include <TopAbs.hxx>
#include <TopoDS_Iterator.hxx>
#include <algorithm>

namespace BRepExt
{
//! Carry intent through Boolean edge modifications, retaining only final edges.
template<class Maker>
static std::vector<blending::Join> rangeJoins(const std::vector<blending::Join>& before,
                                            Maker& maker)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> retained;
  TopExp::MapShapes(maker.Shape(), TopAbs_EDGE, retained);
  std::vector<blending::Join> after;
  for (const auto& join : before)
  {
    if (retained.Contains(join.edge))
      after.push_back(join);
    for (const auto& edge : maker.Modified(join.edge))
      if (edge.ShapeType() == TopAbs_EDGE && retained.Contains(edge))
        after.push_back({edge, join.role, join.order});
  }
  return after;
}

//! True when `shape` is one solid bounded by one shell.
static bool oneShell(const TopoDS_Shape& shape)
{
  int solids = 0, shells = 0;
  for (TopExp_Explorer solid(shape, TopAbs_SOLID); solid.More(); solid.Next())
    ++solids;
  for (TopExp_Explorer shell(shape, TopAbs_SHELL); shell.More(); shell.Next())
    ++shells;
  return solids == 1 && shells == 1;
}

//! Assemble the signed boundary difference from support-face provenance.
//! Unchanged faces cancel. The remaining skins share one intersection inventory;
//! new sheets enter with reversed orientation in the removed-material boundary.
//! Each connected component is then a whole outer shell of the removed or the added
//! material, told apart by its signed volume. That holds when both bodies are single
//! shells sharing an unchanged face: neither difference then has a cavity shell
//! (derived; see the operations documentation).
struct RangeDelta
{
  TopoDS_Compound             removed, added;
  std::vector<blending::Join> joins;
};

static RangeDelta localDelta(const char*              op,
                             const TopoDS_Shape&      shape,
                             const blending::Blended& made)
{
  RangeDelta   delta;
  BRep_Builder builder;
  builder.MakeCompound(delta.removed);
  builder.MakeCompound(delta.added);
  TopoDS_Compound oldSkin, newSkin;
  builder.MakeCompound(oldSkin);
  builder.MakeCompound(newSkin);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> cancelled;
  for (const auto& [old, now] : made.unchanged)
  {
    cancelled.Add(old);
    cancelled.Add(now);
  }
  for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
    if (!cancelled.Contains(face.Current()))
      builder.Add(oldSkin, face.Current());
  for (TopExp_Explorer face(made.shape, TopAbs_FACE); face.More(); face.Next())
    if (!cancelled.Contains(face.Current()))
      builder.Add(newSkin, face.Current());

  // Both oriented boundary differences must share the split edge inventory.
  NCollection_List<TopoDS_Shape> arguments;
  arguments.Append(oldSkin);
  arguments.Append(newSkin);
  BOPAlgo_PaveFiller intersections;
  intersections.SetArguments(arguments);
  intersections.SetNonDestructive(true);
  intersections.Perform();
  BRepAlgoAPI_Cut oldOnly(oldSkin, newSkin, intersections);
  BRepAlgoAPI_Cut newOnly(newSkin, oldSkin, intersections);
  if (!oldOnly.IsDone() || !newOnly.IsDone())
    fail(op, no_solution, "support boundary difference failed");

  // Substitute only retained source-face pieces. No sewing is allowed to
  // rewrite the surfaces or the edge representation lists during assembly.
  BRepTools_ReShape reshape;
  for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
  {
    if (cancelled.Contains(face.Current()) || oldOnly.IsDeleted(face.Current()))
      reshape.Remove(face.Current());
    else if (!oldOnly.Modified(face.Current()).IsEmpty())
    {
      TopoDS_Compound images;
      builder.MakeCompound(images);
      for (const auto& image : oldOnly.Modified(face.Current()))
        if (image.ShapeType() == TopAbs_FACE)
          builder.Add(images, image);
      reshape.Replace(face.Current(), images);
    }
  }
  std::vector<TopoDS_Shape> patches;
  const auto                retained = reshape.Apply(shape);
  for (TopExp_Explorer face(retained, TopAbs_FACE); face.More(); face.Next())
    patches.push_back(face.Current());
  for (TopExp_Explorer face(newOnly.Shape(), TopAbs_FACE); face.More(); face.Next())
    patches.push_back(face.Current().Reversed());

  // Connected components are determined by topological edges, never by a
  // tolerance merge. Each component carries the signed boundary orientation.
  std::vector<bool> used(patches.size(), false);
  for (size_t seed = 0; seed < patches.size(); ++seed)
  {
    if (used[seed])
      continue;
    TopoDS_Shell shell;
    builder.MakeShell(shell);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
    std::vector<size_t>                                        pending{seed};
    used[seed] = true;
    for (size_t cursor = 0; cursor < pending.size(); ++cursor)
    {
      const auto& face = patches[pending[cursor]];
      builder.Add(shell, face);
      TopExp::MapShapes(face, TopAbs_EDGE, edges);
      for (size_t i = 0; i < patches.size(); ++i)
      {
        if (used[i])
          continue;
        for (TopExp_Explorer edge(patches[i], TopAbs_EDGE); edge.More(); edge.Next())
          if (edges.Contains(edge.Current()))
          {
            used[i] = true;
            pending.push_back(i);
            break;
          }
      }
    }
    auto solid = BRepBuilderAPI_MakeSolid(shell).Solid();
    if (!BRepCheck_Analyzer(solid).IsValid())
      fail(op, no_solution, "material boundary is not closed");
    GProp_GProps mass;
    BRepGProp::VolumeProperties(solid, mass);
    const bool removed = mass.Mass() > 0;
    if (!removed)
      solid.Reverse();
    builder.Add(removed ? delta.removed : delta.added, solid);
  }
  delta.joins         = rangeJoins(made.joins, oldOnly);
  const auto newJoins = rangeJoins(made.joins, newOnly);
  delta.joins.insert(delta.joins.end(), newJoins.begin(), newJoins.end());
  return delta;
}

//! A plane square to the edge at a limit of the range: through `point`, `normal` along the range.
struct Limit
{
  gp_Pnt point;
  gp_Dir normal;
  double parameter; //!< the edge curve's parameter there
};

static bool insideSlab(const Limit (&limit)[2], const gp_Pnt& p, double tol)
{
  return gp_Vec(limit[0].point, p).Dot(gp_Vec(limit[0].normal)) > tol
         && gp_Vec(limit[1].point, p).Dot(gp_Vec(limit[1].normal)) < -tol;
}

static double offPlane(const Limit& limit, const gp_Pnt& p)
{
  return gp_Vec(limit.point, p).Dot(gp_Vec(limit.normal));
}

using EdgeFaces =
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>;

//! What the direct limited application is built from (see directPremise).
struct DirectPlan
{
  TopoDS_Face faces[2];    //!< the edge's two faces in the body, as its shell holds them
  TopoDS_Face images[2];   //!< their one image each in the blended body
  TopoDS_Face sheet;       //!< the blend's one new face, as the blended body holds it
  TopoDS_Edge contacts[2]; //!< the sheet's edge on images[i]
  double      crossing[2][2]; //!< contacts[i]'s parameter on limit k
  TopoDS_Shape shell;      //!< the body's shell holding both faces, as the body holds it
};

//! Root of the plane's signed distance along `curve` within [lo, hi], where it changes sign:
//! bisection kept by Newton steps that stay inside the bracket.
static double crossingOf(const Adaptor3d_Curve& curve, const Limit& limit, double lo, double hi)
{
  double glo = offPlane(limit, curve.Value(lo));
  for (int i = 0; i < 200 && hi - lo > 1e-15 * (1 + std::abs(lo) + std::abs(hi)); ++i)
  {
    const double mid = 0.5 * (lo + hi);
    gp_Pnt       p;
    gp_Vec       d;
    curve.D1(mid, p, d);
    const double g = offPlane(limit, p), slope = d.Dot(gp_Vec(limit.normal));
    if (g == 0)
      return mid;
    double next = slope != 0 ? mid - g / slope : mid;
    if ((g < 0) == (glo < 0))
      lo = mid, glo = g;
    else
      hi = mid;
    if (next > lo && next < hi)
    {
      const double gn = offPlane(limit, curve.Value(next));
      if (gn == 0)
        return next;
      if ((gn < 0) == (glo < 0))
        lo = next, glo = gn;
      else
        hi = next;
    }
  }
  return 0.5 * (lo + hi);
}

//! The premise of the direct limited application, with what it is built from:
//! - the blend is this module's own, of one edge with two distinct faces, each changed into one
//!   face, and the blend made one new face (the sheet) meeting each image in one edge (contact);
//! - both limits lie strictly inside the edge;
//! - every edge the blend changed or made other than the contacts lies strictly outside the
//!   slab between the limits, and each contact enters and leaves the slab once, through the two
//!   limits, inside its parameter range.
//! Then within the slab the change is the strip of each face between the edge and its contact
//! and the sheet between the contacts, and each limit meets it in three curves: across each
//! face from the edge to the contact and across the sheet between the contacts. Sampled as the
//! range's own checks are; not an interval certificate.
static bool directPremise(const TopoDS_Shape&      shape,
                          const blending::Blended& made,
                          const TopoDS_Edge&       edge,
                          double                   start,
                          double                   end,
                          const Limit (&limit)[2],
                          DirectPlan&              plan)
{
  if (start <= 0 || end >= 1 || !BRep_Tool::SameParameter(edge))
    return false;
  EdgeFaces bodyFaces;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, bodyFaces);
  const auto* around = bodyFaces.Seek(edge);
  if (!around || around->Extent() != 2 || around->First().IsSame(around->Last()))
    return false;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> kept, imaged;
  for (const auto& [old, now] : made.unchanged)
    kept.Add(now);
  for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
    for (const auto& image : made.history->Modified(face.Current()))
      imaged.Add(image);
  for (int i = 0; i < 2; ++i)
  {
    const auto& images = made.history->Modified(i ? around->Last() : around->First());
    if (images.Extent() != 1 || images.First().ShapeType() != TopAbs_FACE)
      return false;
    plan.images[i] = TopoDS::Face(images.First());
  }
  int sheets = 0;
  for (TopExp_Explorer face(made.shape, TopAbs_FACE); face.More(); face.Next())
    if (!imaged.Contains(face.Current()) && !kept.Contains(face.Current()))
    {
      plan.sheet = TopoDS::Face(face.Current());
      ++sheets;
    }
  if (sheets != 1)
    return false;
  EdgeFaces madeFaces;
  TopExp::MapShapesAndAncestors(made.shape, TopAbs_EDGE, TopAbs_FACE, madeFaces);
  int found[2] = {0, 0};
  for (TopExp_Explorer e(plan.sheet, TopAbs_EDGE); e.More(); e.Next())
    for (int i = 0; i < 2; ++i)
      for (const auto& face : madeFaces.FindFromKey(e.Current()))
        if (face.IsSame(plan.images[i]) && !plan.contacts[i].IsSame(e.Current()))
        {
          plan.contacts[i] = TopoDS::Edge(e.Current());
          ++found[i];
        }
  if (found[0] != 1 || found[1] != 1 || plan.contacts[0].IsSame(plan.contacts[1]))
    return false;
  for (const auto& contact : plan.contacts)
    if (BRep_Tool::IsClosed(contact, plan.sheet) || !BRep_Tool::SameParameter(contact))
      return false;
  // Everything else the blend changed or made stays out of the slab.
  auto outside = [&](const TopoDS_Shape& candidate) {
    BRepAdaptor_Curve curve(TopoDS::Edge(candidate));
    for (int i = 0; i <= 64; ++i)
      if (insideSlab(limit,
                     curve.Value(curve.FirstParameter()
                                 + (curve.LastParameter() - curve.FirstParameter()) * i / 64),
                     -BRepExt::tolerance::length))
        return false;
    return true;
  };
  for (const TopoDS_Face& face : {plan.sheet, plan.images[0], plan.images[1]})
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
      if (!kept.Contains(e.Current()) && !e.Current().IsSame(plan.contacts[0])
          && !e.Current().IsSame(plan.contacts[1]) && !BRep_Tool::Degenerated(TopoDS::Edge(e.Current()))
          && !outside(e.Current()))
        return false;
  // Each contact runs through the slab once, the way the edge does.
  BRepAdaptor_Curve along(edge);
  const double      middle = along.FirstParameter() + (along.LastParameter() - along.FirstParameter())
                                                   * (edge.Orientation() == TopAbs_REVERSED ? 1 - (start + end) / 2
                                                                                            : (start + end) / 2);
  gp_Pnt mid;
  gp_Vec forward;
  along.D1(middle, mid, forward);
  if (limit[1].parameter < limit[0].parameter)
    forward.Reverse();
  for (int i = 0; i < 2; ++i)
  {
    BRepAdaptor_Curve contact(plan.contacts[i]);
    const double      f = contact.FirstParameter(), l = contact.LastParameter();
    constexpr int     samples = 256;
    int               runs = 0, first = -1, last = -1;
    bool              was = false;
    for (int k = 0; k <= samples; ++k)
    {
      const bool in = insideSlab(limit, contact.Value(f + (l - f) * k / samples), 0);
      if (in && !was)
        ++runs, first = k;
      if (in)
        last = k;
      was = in;
    }
    if (runs != 1 || first == 0 || last == samples)
      return false;
    gp_Pnt p;
    gp_Vec d;
    contact.D1(f + (l - f) * (first + last) / 2 / samples, p, d);
    const bool codirected = d.Dot(forward) > 0;
    const int  entry = codirected ? 0 : 1;
    const double into[2] = {f + (l - f) * (first - 1) / samples, f + (l - f) * first / samples};
    const double outof[2] = {f + (l - f) * last / samples, f + (l - f) * (last + 1) / samples};
    for (const auto& [k, bracket] : {std::pair{entry, into}, std::pair{1 - entry, outof}})
    {
      const double g0 = offPlane(limit[k], contact.Value(bracket[0])),
                   g1 = offPlane(limit[k], contact.Value(bracket[1]));
      if ((g0 < 0) == (g1 < 0) && g0 != 0 && g1 != 0)
        return false;
      plan.crossing[i][k] = crossingOf(contact, limit[k], bracket[0], bracket[1]);
    }
  }
  // The faces as the body's shell holds them (an unplaced shell: faces are put back in it).
  for (TopExp_Explorer shell(shape, TopAbs_SHELL); shell.More() && plan.shell.IsNull(); shell.Next())
  {
    int held = 0;
    for (TopoDS_Iterator face(shell.Current().Oriented(TopAbs_FORWARD)); face.More(); face.Next())
      for (int i = 0; i < 2; ++i)
        if (face.Value().IsSame(i ? around->Last() : around->First()))
        {
          plan.faces[i] = TopoDS::Face(face.Value());
          ++held;
        }
    if (held == 2)
      plan.shell = shell.Current();
  }
  return !plan.shell.IsNull() && plan.shell.Location().IsIdentity();
}

//! One Newton solve across: the step `lambda` along `across` from `q` in the surface's
//! parameters to the plane.
static bool solveAcross(const Adaptor3d_Surface& surface,
                        const Limit&             limit,
                        const gp_Pnt2d&          q,
                        const gp_Vec2d&          across,
                        double&                  lambda)
{
  for (int i = 0; i < 50; ++i)
  {
    gp_Pnt p;
    gp_Vec du, dv;
    surface.D1(q.X() + lambda * across.X(), q.Y() + lambda * across.Y(), p, du, dv);
    const double g = offPlane(limit, p), slope = (du * across.X() + dv * across.Y()).Dot(gp_Vec(limit.normal));
    if (std::abs(g) <= BRepExt::tolerance::rounding(gp_Vec(p.XYZ()).Magnitude()))
      return true;
    if (slope == 0)
      return false;
    lambda -= g / slope;
  }
  return false;
}

//! A curve in a plane, in that plane's parameters: lines and conics by projection, which is exact
//! for them; splines by their poles, which an affine map carries exactly.
static occ::handle<Geom2d_Curve> inPlane(const occ::handle<Geom_Curve>& curve,
                                        double                         first,
                                        double                         last,
                                        const occ::handle<Geom_Plane>& plane)
{
  if (const auto spline = occ::down_cast<Geom_BSplineCurve>(curve))
  {
    NCollection_Array1<gp_Pnt2d> poles(1, spline->NbPoles());
    for (int i = 1; i <= spline->NbPoles(); ++i)
    {
      double u, v;
      ElSLib::Parameters(plane->Pln(), spline->Pole(i), u, v);
      poles(i) = gp_Pnt2d(u, v);
    }
    if (spline->IsRational())
      return new Geom2d_BSplineCurve(poles, *spline->Weights(), spline->Knots(), spline->Multiplicities(),
                                     spline->Degree(), spline->IsPeriodic());
    return new Geom2d_BSplineCurve(poles, spline->Knots(), spline->Multiplicities(), spline->Degree(),
                                   spline->IsPeriodic());
  }
  auto projected = GeomProjLib::Curve2d(curve, first, last, plane);
  if (const auto trimmed = occ::down_cast<Geom2d_TrimmedCurve>(projected))
    projected = trimmed->BasisCurve();
  return projected;
}

//! The curve where a limit's plane crosses `face` from the point `from` (parameters `a`) to `to`
//! (parameters `b`), its curve on the face and on `plane`, running from `from` to `to`.
//! Elementary surfaces meet a plane in lines and conics, taken exact (OpenCascade's analytic
//! intersection); other surfaces by one solve across per sample (solveAcross), an isoparametric
//! curve where it lies in the plane to rounding, else interpolated in both parameter spaces at the
//! same parameters, refined while the curve on the surface and in the plane come closer.
struct Section
{
  occ::handle<Geom_Curve>   curve;
  occ::handle<Geom2d_Curve> onFace, onPlane;
  double                    first = 0, last = 0, gap = 0;
};

static Section planeSection(const char*                    op,
                            const TopoDS_Face&             face,
                            const Limit&                   limit,
                            const occ::handle<Geom_Plane>& plane,
                            const gp_Pnt2d&                a,
                            const gp_Pnt2d&                b,
                            const gp_Pnt&                  from,
                            const gp_Pnt&                  to)
{
  const auto          surface = BRep_Tool::Surface(face);
  GeomAdaptor_Surface adaptor(surface);
  gp_Vec2d            across(-(b.Y() - a.Y()), b.X() - a.X());
  if (across.Magnitude() == 0)
    fail(op, no_solution, "a limit's section has no extent");
  across.Normalize();
  Section section;
  auto measure = [&](int samples) {
    double gap = 0;
    for (int k = 0; k <= samples; ++k)
    {
      const double t  = section.first + (section.last - section.first) * k / samples;
      const auto   p  = section.curve->Value(t);
      const auto   uv = section.onFace->Value(t);
      const auto   w  = section.onPlane->Value(t);
      gap = std::max({gap, p.Distance(adaptor.Value(uv.X(), uv.Y())), p.Distance(plane->Value(w.X(), w.Y()))});
    }
    return gap;
  };
  const auto type = adaptor.GetType();
  if (type == GeomAbs_Plane || type == GeomAbs_Cylinder || type == GeomAbs_Cone || type == GeomAbs_Sphere
      || type == GeomAbs_Torus)
  {
    GeomAPI_IntSS crossing(plane, surface, BRepExt::tolerance::length);
    if (!crossing.IsDone())
      fail(op, no_solution, "a limit's plane does not cross a face");
    double best = RealLast();
    for (int i = 1; i <= crossing.NbLines(); ++i)
    {
      auto line = crossing.Line(i);
      if (const auto trimmed = occ::down_cast<Geom_TrimmedCurve>(line))
        line = trimmed->BasisCurve();
      GeomAPI_ProjectPointOnCurve onA(from, line), onB(to, line);
      if (!onA.NbPoints() || !onB.NbPoints())
        continue;
      const double off = std::max(onA.LowerDistance(), onB.LowerDistance());
      if (off < best)
      {
        best           = off;
        section.curve  = line;
        section.first  = onA.LowerDistanceParameter();
        section.last   = onB.LowerDistanceParameter();
      }
    }
    if (section.curve.IsNull())
      fail(op, no_solution, "a limit's plane does not cross a face");
    if (section.curve->IsPeriodic())
    {
      // The arc through the section's middle, found by the solve across from the middle of
      // the parameter segment.
      const double period = section.curve->Period();
      double       lambda = 0;
      const gp_Pnt2d q((a.X() + b.X()) / 2, (a.Y() + b.Y()) / 2);
      if (!solveAcross(adaptor, limit, q, across, lambda))
        fail(op, no_solution, "a limit's section does not converge");
      GeomAPI_ProjectPointOnCurve onMiddle(adaptor.Value(q.X() + lambda * across.X(), q.Y() + lambda * across.Y()),
                                           section.curve);
      double middle = onMiddle.LowerDistanceParameter();
      auto   above  = [&](double t, double base) {
        while (t < base)
          t += period;
        while (t >= base + period)
          t -= period;
        return t;
      };
      double last = above(section.last, section.first);
      if (above(middle, section.first) > last)
      {
        section.curve = section.curve->Reversed();
        section.first = section.curve->ReversedParameter(section.first);
        last          = above(section.curve->ReversedParameter(section.last), section.first);
      }
      section.last = last;
    }
    else if (section.last < section.first)
    {
      section.curve = section.curve->Reversed();
      section.first = section.curve->ReversedParameter(section.first);
      section.last  = section.curve->ReversedParameter(section.last);
    }
    section.onFace = GeomProjLib::Curve2d(section.curve, section.first, section.last, surface);
    if (const auto trimmed = occ::down_cast<Geom2d_TrimmedCurve>(section.onFace))
      section.onFace = trimmed->BasisCurve();
    section.onPlane = inPlane(section.curve, section.first, section.last, plane);
    if (section.onFace.IsNull() || section.onPlane.IsNull())
      fail(op, no_solution, "a limit's section has no curve on a face");
    section.gap = measure(16);
    return section;
  }
  // Freeform: one solve across per sample.
  double previous = RealLast();
  for (int samples = 16;; samples *= 2)
  {
    NCollection_Array1<gp_Pnt2d> uv(1, samples + 1), w(1, samples + 1);
    auto                         params = new NCollection_HArray1<double>(1, samples + 1);
    double                       lambda = 0;
    for (int k = 0; k <= samples; ++k)
    {
      params->SetValue(k + 1, double(k) / samples);
      if (k == 0 || k == samples)
      {
        uv(k + 1) = k ? b : a;
        double u, v;
        ElSLib::Parameters(plane->Pln(), k ? to : from, u, v);
        w(k + 1) = gp_Pnt2d(u, v);
        continue;
      }
      const gp_Pnt2d q(a.X() + (b.X() - a.X()) * k / samples, a.Y() + (b.Y() - a.Y()) * k / samples);
      if (!solveAcross(adaptor, limit, q, across, lambda))
        fail(op, no_solution, "a limit's section does not converge");
      uv(k + 1) = gp_Pnt2d(q.X() + lambda * across.X(), q.Y() + lambda * across.Y());
      double u, v;
      ElSLib::Parameters(plane->Pln(), adaptor.Value(uv(k + 1).X(), uv(k + 1).Y()), u, v);
      w(k + 1) = gp_Pnt2d(u, v);
    }
    // An isoparametric section, where the solved points lie on the surface's iso curve to
    // rounding (the ends are the contacts' points, which their tolerance covers).
    for (const bool alongV : {false, true})
    {
      const double iso = alongV ? (a.X() + b.X()) / 2 : (a.Y() + b.Y()) / 2;
      bool         flat = true;
      for (int k = 2; k <= samples && flat; ++k)
      {
        const auto p = alongV ? adaptor.Value(iso, uv(k).Y()) : adaptor.Value(uv(k).X(), iso);
        flat = p.Distance(adaptor.Value(uv(k).X(), uv(k).Y()))
               <= BRepExt::tolerance::rounding(gp_Vec(p.XYZ()).Magnitude());
      }
      if (!flat)
        continue;
      section.curve = alongV ? surface->UIso(iso) : surface->VIso(iso);
      section.first = alongV ? a.Y() : a.X();
      section.last  = alongV ? b.Y() : b.X();
      // On the face the line t -> (iso, t); reversed, the curve's parameter is r(t) = c - t, so
      // the line through (iso, r(0)) backwards.
      double sense = 1, origin = 0;
      if (section.last < section.first)
      {
        origin         = section.curve->ReversedParameter(0);
        sense          = -1;
        section.curve  = section.curve->Reversed();
        section.first  = section.curve->ReversedParameter(section.first);
        section.last   = section.curve->ReversedParameter(section.last);
      }
      section.onFace = alongV ? new Geom2d_Line(gp_Pnt2d(iso, origin), gp_Dir2d(0, sense))
                              : new Geom2d_Line(gp_Pnt2d(origin, iso), gp_Dir2d(sense, 0));
      section.onPlane = inPlane(section.curve, section.first, section.last, plane);
      section.gap     = measure(2 * samples);
      return section;
    }
    auto pointsOn = new NCollection_HArray1<gp_Pnt2d>(uv), pointsIn = new NCollection_HArray1<gp_Pnt2d>(w);
    Geom2dAPI_Interpolate onFace(pointsOn, params, false, Precision::PConfusion()),
      onPlane(pointsIn, params, false, Precision::PConfusion());
    onFace.Perform();
    onPlane.Perform();
    if (!onFace.IsDone() || !onPlane.IsDone())
      fail(op, no_solution, "a limit's section cannot be interpolated");
    Section trial;
    trial.onFace  = onFace.Curve();
    trial.onPlane = onPlane.Curve();
    trial.curve   = GeomAPI::To3d(trial.onPlane, plane->Pln());
    trial.first   = 0;
    trial.last    = 1;
    std::swap(section, trial);
    section.gap = measure(4 * samples);
    if (section.gap <= BRepExt::tolerance::length)
      return section;
    if (section.gap >= previous)
    {
      std::swap(section, trial); // the finer fit came no closer: keep the coarser
      return section;
    }
    previous = section.gap;
  }
}

//! The direct limited application: within the slab, the edge's faces lose the strips between
//! the edge and the contacts, the sheet between the limits enters, and each limit closes the
//! change with a planar cap; the body's other faces, edges and vertices are kept as they are.
static TopoDS_Shape directRange(const char*              op,
                                const TopoDS_Shape&      shape,
                                const blending::Blended& made,
                                const TopoDS_Edge&       edge,
                                const Limit (&limit)[2],
                                const DirectPlan&        plan,
                                FaceHistory&             history,
                                std::vector<blending::Join>* joins)
{
  BRep_Builder builder;
  const double length = BRepExt::tolerance::length;
  auto vertex = [&](const gp_Pnt& p) {
    TopoDS_Vertex v;
    builder.MakeVertex(v, p, length);
    return v;
  };
  auto bounded = [&](TopoDS_Edge e, const TopoDS_Vertex& v0, const TopoDS_Vertex& v1, double f, double l) {
    builder.Add(e, v0.Oriented(TopAbs_FORWARD));
    builder.Add(e, v1.Oriented(TopAbs_REVERSED));
    builder.Range(e, f, l);
    return e;
  };
  // The edge's pieces outside the slab, on the faces' own curves.
  const auto   forwardEdge = TopoDS::Edge(edge.Oriented(TopAbs_FORWARD));
  double       ef, el;
  BRep_Tool::Range(forwardEdge, ef, el);
  TopoDS_Vertex ev0, ev1;
  TopExp::Vertices(forwardEdge, ev0, ev1);
  const int    low = limit[0].parameter < limit[1].parameter ? 0 : 1; // the limit nearer the first end
  TopoDS_Vertex onEdge[2] = {vertex(limit[0].point), vertex(limit[1].point)};
  const auto   before = bounded(TopoDS::Edge(forwardEdge.EmptyCopied()), ev0, onEdge[low], ef, limit[low].parameter);
  const auto   after  = bounded(TopoDS::Edge(forwardEdge.EmptyCopied()), onEdge[1 - low], ev1, limit[1 - low].parameter, el);
  // The contacts between the limits, on their curves from the blend.
  TopoDS_Vertex onContact[2][2];
  TopoDS_Edge   contact[2];
  for (int i = 0; i < 2; ++i)
  {
    double     f, l;
    TopLoc_Location location;
    const auto curve = BRep_Tool::Curve(plan.contacts[i], location, f, l);
    BRepAdaptor_Curve along(plan.contacts[i]);
    for (int k = 0; k < 2; ++k)
      onContact[i][k] = vertex(along.Value(plan.crossing[i][k]));
    const int first = plan.crossing[i][0] < plan.crossing[i][1] ? 0 : 1;
    builder.MakeEdge(contact[i], curve, location, BRep_Tool::Tolerance(plan.contacts[i]));
    builder.UpdateEdge(contact[i], BRep_Tool::CurveOnSurface(plan.contacts[i], plan.sheet, f, l), plan.sheet,
                       BRep_Tool::Tolerance(plan.contacts[i]));
    builder.UpdateEdge(contact[i], BRep_Tool::CurveOnSurface(plan.contacts[i], plan.images[i], f, l), plan.faces[i],
                       BRep_Tool::Tolerance(plan.contacts[i]));
    const auto regularity = BRep_Tool::Continuity(plan.contacts[i], plan.sheet, plan.images[i]);
    if (regularity != GeomAbs_C0)
      builder.Continuity(contact[i], plan.sheet, plan.faces[i], regularity);
    bounded(contact[i], onContact[i][first], onContact[i][1 - first], plan.crossing[i][first], plan.crossing[i][1 - first]);
  }
  // The caps' planes and the sections across the faces and the sheet.
  occ::handle<Geom_Plane> planes[2] = {new Geom_Plane(gp_Pln(limit[0].point, limit[0].normal)),
                                       new Geom_Plane(gp_Pln(limit[1].point, limit[1].normal))};
  TopoDS_Face caps[2];
  for (int k = 0; k < 2; ++k)
    builder.MakeFace(caps[k], planes[k], length);
  auto sectionEdge = [&](const Section& s, const TopoDS_Face& face, const TopoDS_Face& cap, const TopoDS_Vertex& v0,
                         const TopoDS_Vertex& v1) {
    TopoDS_Edge e;
    const double tol = std::max(length, s.gap);
    builder.MakeEdge(e, s.curve, tol);
    builder.UpdateEdge(e, s.onFace, face, tol);
    builder.UpdateEdge(e, s.onPlane, cap, tol);
    return bounded(e, v0, v1, s.first, s.last);
  };
  TopoDS_Edge across[2][2], sections[2];
  for (int k = 0; k < 2; ++k)
  {
    for (int i = 0; i < 2; ++i)
    {
      double f, l;
      const auto onFace = BRep_Tool::CurveOnSurface(forwardEdge, plan.faces[i], f, l)->Value(limit[k].parameter);
      const auto onImage = BRep_Tool::CurveOnSurface(plan.contacts[i], plan.images[i], f, l)->Value(plan.crossing[i][k]);
      across[k][i] = sectionEdge(planeSection(op, plan.faces[i], limit[k], planes[k], onFace, onImage,
                                              limit[k].point, BRep_Tool::Pnt(onContact[i][k])),
                                 plan.faces[i], caps[k], onEdge[k], onContact[i][k]);
    }
    double f, l;
    const auto a = BRep_Tool::CurveOnSurface(plan.contacts[0], plan.sheet, f, l)->Value(plan.crossing[0][k]);
    const auto b = BRep_Tool::CurveOnSurface(plan.contacts[1], plan.sheet, f, l)->Value(plan.crossing[1][k]);
    sections[k] = sectionEdge(planeSection(op, plan.sheet, limit[k], planes[k], a, b, BRep_Tool::Pnt(onContact[0][k]),
                                           BRep_Tool::Pnt(onContact[1][k])),
                              plan.sheet, caps[k], onContact[0][k], onContact[1][k]);
  }
  // Orient a loop: from `first` as given, each next edge starts where the last one ended.
  auto loop = [&](std::vector<TopoDS_Edge> edges, const TopoDS_Edge& first, bool closed = true) {
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    builder.Add(wire, first);
    auto at = TopExp::LastVertex(first, true);
    edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const TopoDS_Edge& e) { return e.IsSame(first); }),
                edges.end());
    while (!edges.empty())
    {
      auto next = std::find_if(edges.begin(), edges.end(), [&](const TopoDS_Edge& e) {
        return TopExp::FirstVertex(e).IsSame(at) || TopExp::LastVertex(e).IsSame(at);
      });
      if (next == edges.end())
        fail(op, no_solution, "the limited change does not close");
      const auto base     = TopoDS::Edge(next->Oriented(TopAbs_FORWARD));
      const auto oriented = TopExp::FirstVertex(base).IsSame(at) ? base : TopoDS::Edge(base.Reversed());
      builder.Add(wire, oriented);
      at = TopExp::LastVertex(oriented, true);
      edges.erase(next);
    }
    if (closed && !at.IsSame(TopExp::FirstVertex(first, true)))
      fail(op, no_solution, "the limited change does not close");
    wire.Closed(closed);
    return wire;
  };
  auto orientationIn = [](const TopoDS_Shape& face, const TopoDS_Shape& e) {
    for (TopExp_Explorer x(face.Oriented(TopAbs_FORWARD), TopAbs_EDGE); x.More(); x.Next())
      if (x.Current().IsSame(e))
        return x.Current().Orientation();
    return TopAbs_EXTERNAL;
  };
  // The edge's faces: the edge's piece between the limits gives way to the strip's other sides.
  TopoDS_Face changed[2];
  for (int i = 0; i < 2; ++i)
  {
    const auto contactSide = orientationIn(plan.images[i], plan.contacts[i]);
    changed[i] = TopoDS::Face(plan.faces[i].Oriented(TopAbs_FORWARD).EmptyCopied());
    for (TopoDS_Iterator w(plan.faces[i].Oriented(TopAbs_FORWARD)); w.More(); w.Next())
    {
      bool holds = false;
      for (TopoDS_Iterator e(w.Value()); e.More(); e.Next())
        holds |= e.Value().IsSame(edge);
      if (!holds)
      {
        builder.Add(changed[i], w.Value());
        continue;
      }
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      for (TopoDS_Iterator e(w.Value()); e.More(); e.Next())
      {
        if (!e.Value().IsSame(edge))
        {
          builder.Add(wire, e.Value());
          continue;
        }
        const bool forward = e.Value().Orientation() == TopAbs_FORWARD;
        // Forward: first end, low limit, contact, high limit, last end.
        const auto pieceLow = forward ? before : TopoDS::Edge(after.Reversed());
        auto       path     = loop({pieceLow, across[forward ? low : 1 - low][i], contact[i],
                                    across[forward ? 1 - low : low][i], forward ? after : TopoDS::Edge(before.Reversed())},
                                   pieceLow, false);
        for (TopoDS_Iterator p(path); p.More(); p.Next())
        {
          if (p.Value().IsSame(contact[i]) && p.Value().Orientation() != contactSide)
            fail(op, no_solution, "the contact runs against its face");
          builder.Add(wire, p.Value());
        }
      }
      // The iterator composed the wire's orientation into its edges.
      wire.Closed(true);
      builder.Add(changed[i], wire);
    }
    changed[i].Orientation(plan.faces[i].Orientation());
  }
  // The sheet between the limits, as the blend oriented it.
  TopAbs_Orientation sheetSide = TopAbs_FORWARD;
  for (TopExp_Explorer f(made.shape, TopAbs_FACE); f.More(); f.Next())
    if (f.Current().IsSame(plan.sheet))
      sheetSide = f.Current().Orientation();
  auto piece = TopoDS::Face(plan.sheet.Oriented(TopAbs_FORWARD).EmptyCopied());
  {
    const auto first = orientationIn(plan.sheet, plan.contacts[0]) == TopAbs_FORWARD ? contact[0]
                                                                                    : TopoDS::Edge(contact[0].Reversed());
    const auto wire = loop({first, sections[0], contact[1], sections[1]}, first);
    for (TopoDS_Iterator p(wire); p.More(); p.Next())
      if (p.Value().IsSame(contact[1]) && p.Value().Orientation() != orientationIn(plan.sheet, plan.contacts[1]))
        fail(op, no_solution, "the sheet's contacts run against each other");
    builder.Add(piece, wire);
  }
  // The body's shell orientation: faces are placed relative to it.
  const auto shellSide = plan.shell.Orientation();
  auto relative = [&](TopAbs_Orientation absolute) {
    return shellSide == TopAbs_REVERSED ? TopAbs::Reverse(absolute) : absolute;
  };
  piece.Orientation(relative(sheetSide));
  // Each cap: the section used against the sheet's use of it; the face forward where its loop
  // turns positively round the plane's normal.
  for (int k = 0; k < 2; ++k)
  {
    const auto sheetUse = orientationIn(piece.Oriented(TopAbs_FORWARD), sections[k]);
    const auto absoluteSheetUse = sheetSide == TopAbs_REVERSED ? TopAbs::Reverse(sheetUse) : sheetUse;
    const auto first = absoluteSheetUse == TopAbs_FORWARD ? TopoDS::Edge(sections[k].Reversed()) : sections[k];
    auto wire = loop({first, across[k][0], across[k][1]}, first);
    double area = 0;
    gp_Pnt2d previous;
    bool started = false;
    for (TopoDS_Iterator p(wire); p.More(); p.Next())
    {
      const auto e = TopoDS::Edge(p.Value());
      double     f, l;
      const auto pc = BRep_Tool::CurveOnSurface(e, caps[k], f, l);
      for (int s = 0; s <= 16; ++s)
      {
        const double t = e.Orientation() == TopAbs_FORWARD ? f + (l - f) * s / 16 : l - (l - f) * s / 16;
        const auto   w = pc->Value(t);
        if (started)
          area += previous.X() * w.Y() - w.X() * previous.Y();
        previous = w;
        started  = true;
      }
    }
    // Absolute orientation of the cap: forward where the loop turns positively.
    builder.Add(caps[k], area > 0 ? wire : TopoDS::Wire(wire.Reversed()));
    caps[k].Orientation(relative(area > 0 ? TopAbs_FORWARD : TopAbs_REVERSED));
  }
  // Tolerances from the real gaps: each new vertex covers the ends of its edges' curves.
  double allowed = length;
  for (const auto& e : {contact[0], contact[1], sections[0], sections[1], across[0][0], across[0][1],
                        across[1][0], across[1][1], before, after})
    allowed = std::max(allowed, BRep_Tool::Tolerance(e));
  auto cover = [&](const TopoDS_Edge& e, const std::vector<TopoDS_Face>& faces) {
    double f, l;
    TopoDS_Vertex v[2];
    TopExp::Vertices(e, v[0], v[1]);
    const auto curve = BRep_Tool::Curve(e, f, l);
    for (int end = 0; end < 2; ++end)
    {
      const auto p   = BRep_Tool::Pnt(v[end]);
      double     gap = std::max(BRep_Tool::Tolerance(e), p.Distance(curve->Value(end ? l : f)));
      for (const auto& face : faces)
      {
        double     pf, pl;
        const auto pc = BRep_Tool::CurveOnSurface(e, face, pf, pl);
        const auto uv = pc->Value(end ? pl : pf);
        gap = std::max(gap, p.Distance(BRep_Tool::Surface(face)->Value(uv.X(), uv.Y())));
      }
      // The ends meet within the tolerances the curves already carry (the contacts' and the
      // sections'); a wider gap is a construction that does not close.
      if (gap > allowed)
        fail(op, no_solution, "the limited change's curves do not meet");
      if (gap > BRep_Tool::Tolerance(v[end]))
        builder.UpdateVertex(v[end], gap);
    }
  };
  for (int i = 0; i < 2; ++i)
  {
    cover(contact[i], {piece, changed[i]});
    for (int k = 0; k < 2; ++k)
      cover(across[k][i], {changed[i], caps[k]});
  }
  for (int k = 0; k < 2; ++k)
    cover(sections[k], {piece, caps[k]});
  cover(before, {changed[0], changed[1]});
  cover(after, {changed[0], changed[1]});
  // Assembly by local replacement in the body's shell.
  TopoDS_Shell shell;
  builder.MakeShell(shell);
  for (TopoDS_Iterator f(plan.shell.Oriented(TopAbs_FORWARD)); f.More(); f.Next())
  {
    if (f.Value().IsSame(plan.faces[0]))
      builder.Add(shell, changed[0]);
    else if (f.Value().IsSame(plan.faces[1]))
      builder.Add(shell, changed[1]);
    else
      builder.Add(shell, f.Value());
  }
  builder.Add(shell, piece);
  builder.Add(shell, caps[0]);
  builder.Add(shell, caps[1]);
  shell.Closed(plan.shell.Closed());
  BRepTools_ReShape reshape;
  reshape.Replace(plan.shell.Oriented(TopAbs_FORWARD), shell);
  const auto result = reshape.Apply(shape);
  // History: the edge's faces become their changed copies, the others stay.
  for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
  {
    if (!history.tracks(f.Current()))
      continue;
    NCollection_List<TopoDS_Shape> successors;
    if (f.Current().IsSame(plan.faces[0]))
      successors.Append(changed[0]);
    else if (f.Current().IsSame(plan.faces[1]))
      successors.Append(changed[1]);
    else
      successors.Append(f.Current());
    history.modified(f.Current(), successors);
  }
  if (joins)
  {
    joins->clear();
    for (int i = 0; i < 2; ++i)
      for (const auto& join : made.joins)
        if (join.edge.IsSame(plan.contacts[i]))
        {
          joins->push_back({contact[i], join.role, join.order});
          break;
        }
    for (int k = 0; k < 2; ++k)
      for (const auto& e : {sections[k], across[k][0], across[k][1]})
        joins->push_back({e, blending::JoinRole::Stop, 0});
  }
  return result;
}

//! The part of the change `made` makes to `shape` that lies between the planes square to `edge`
//! at the fractions `start` and `end` of it (the contour's other edges, `others`, must not come
//! back between them): the material removed and added there only, the rest of `shape` as it was.
//! The faces' history goes to `history`; the result is stored.
static EvolutionData withinRange(OcctKernel&                     kernel,
                                 const char*                     op,
                                 const TopoDS_Shape&             shape,
                                 const TopoDS_Shape&             made,
                                 const TopoDS_Edge&              edge,
                                 double                          start,
                                 double                          end,
                                 const std::vector<TopoDS_Edge>& others,
                                 FaceHistory&                    history,
                                 std::vector<blending::Join>*     joins = nullptr,
                                 const blending::Blended*         provenance = nullptr)
{
  BRepAdaptor_Curve curve(edge);
  double a = curve.FirstParameter(), b = curve.LastParameter();
  if (edge.Orientation() == TopAbs_REVERSED)
    std::swap(a, b);
  gp_Pnt p0, p1;
  gp_Vec t0, t1;
  curve.D1(a + (b - a) * start, p0, t0);
  curve.D1(a + (b - a) * end, p1, t1);
  if (t0.Magnitude() < 1e-12 || t1.Magnitude() < 1e-12)
    fail(op, unsupported_geometry, "a limit has a singular tangent");
  if (b < a)
  {
    t0.Reverse();
    t1.Reverse();
  }
  const gp_Dir n0(t0), n1(t1);
  auto         inside = [&](const gp_Pnt& p, double tol) {
    return gp_Vec(p0, p).Dot(gp_Vec(n0)) > tol && gp_Vec(p1, p).Dot(gp_Vec(n1)) < -tol;
  };
  // A pair of normal stopping planes cannot describe a winding interval.
  // Refuse those inputs instead of accidentally rounding a second arc.
  for (int i = 0; i <= 256; ++i)
  {
    const double f = double(i) / 256;
    const auto   p = curve.Value(a + (b - a) * f);
    if (f > start + 1e-8 && f < end - 1e-8)
    {
      if (!inside(p, -BRepExt::tolerance::length))
        fail(op, unsupported_geometry, "range crosses its normal stopping planes");
    }
    else if ((f < start - 1e-8 || f > end + 1e-8) && inside(p, BRepExt::tolerance::length))
      fail(op, unsupported_geometry, "range returns inside its normal stopping planes");
  }
  for (const TopoDS_Edge& onward : others)
  {
    BRepAdaptor_Curve other(onward);
    for (int i = 1; i < 64; ++i)
      if (inside(other.Value(other.FirstParameter()
                             + (other.LastParameter() - other.FirstParameter()) * i / 64),
                 BRepExt::tolerance::length))
        fail(op, unsupported_geometry, "tangent continuation re-enters the limited interval");
  }
  // Where the premise holds, build the limited application directly.
  if (provenance)
  {
    const Limit limit[2] = {{p0, n0, a + (b - a) * start}, {p1, n1, a + (b - a) * end}};
    DirectPlan  plan;
    if (directPremise(shape, *provenance, edge, start, end, limit, plan))
    {
      const auto result = directRange(op, shape, *provenance, edge, limit, plan, history, joins);
      if (!BRepExt::invalidity(result, shape).empty())
        fail(op, no_solution, "limited fillet is invalid");
      return history.finish(store(kernel, result));
    }
  }
  // Isolate both sides of the change. Concave blends add material; convex
  // blends remove it. Mixed contours can do both in one operation.
  // Use the explicit boundary provenance where the blend owns it and the
  // signed-shell ownership is proven. Otherwise, and for the stock builder,
  // the two differences share one safe intersection data structure.
  TopoDS_Shape removedShape, addedShape;
  std::vector<blending::Join> removedIntent, addedIntent;
  if (provenance && !provenance->unchanged.empty() && oneShell(shape) && oneShell(made))
  {
    const auto delta = localDelta(op, shape, *provenance);
    removedShape = delta.removed;
    addedShape = delta.added;
    removedIntent = delta.joins;
    addedIntent = delta.joins;
  }
  else
  {
    NCollection_List<TopoDS_Shape> arguments;
    arguments.Append(shape);
    arguments.Append(made);
    BOPAlgo_PaveFiller intersections;
    intersections.SetArguments(arguments);
    intersections.SetNonDestructive(true);
    intersections.Perform();
    BRepAlgoAPI_Cut removed(shape, made, intersections);
    BRepAlgoAPI_Cut added(made, shape, intersections);
    if (!removed.IsDone() || !added.IsDone())
      fail(op, no_solution, "could not isolate fillet volume");
    removedShape = removed.Shape();
    addedShape = added.Shape();
    if (joins)
    {
      removedIntent = rangeJoins(*joins, removed);
      addedIntent = rangeJoins(*joins, added);
    }
  }
  auto clip = [&](TopoDS_Shape tool, std::vector<blending::Join>& intent) {
    for (int i = 0; i < 2; ++i)
    {
      const auto         point   = i ? p1 : p0;
      const auto         normal  = i ? n1 : n0;
      const auto         plane   = BRepBuilderAPI_MakeFace(gp_Pln(point, normal)).Face();
      const auto         witness = point.Translated(gp_Vec(normal) * (i ? -1 : 1));
      const auto         half    = BRepPrimAPI_MakeHalfSpace(plane, witness).Solid();
      BRepAlgoAPI_Common common(tool, half);
      if (!common.IsDone())
        fail(op, no_solution, "fillet limit intersection failed");
      if (joins)
      {
        intent = rangeJoins(intent, common);
        // The exact plane supplied to Common owns these new faces. They close
        // the limited material change sharply; they are not support contacts.
        const auto support = BRep_Tool::Surface(plane);
        for (TopExp_Explorer f(common.Shape(), TopAbs_FACE); f.More(); f.Next())
          if (BRep_Tool::Surface(TopoDS::Face(f.Current())) == support)
            for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
              intent.push_back({e.Current(), blending::JoinRole::Stop, 0});
      }
      tool = common.Shape();
    }
    return tool;
  };
  Untouched<BRepAlgoAPI_Cut> cut(shape, clip(removedShape, removedIntent));
  if (!cut.IsDone())
    fail(op, no_solution, "limited fillet subtraction failed");
  BRepAlgoAPI_Fuse fuse(cut.Shape(), clip(addedShape, addedIntent));
  if (!fuse.IsDone() || !BRepExt::invalidity(fuse.Shape(), shape).empty())
    fail(op, no_solution, "limited fillet is invalid");
  history.collectComposed(cut, fuse, shape);
  if (joins)
  {
    *joins = rangeJoins(removedIntent, cut);
    joins->insert(joins->end(), addedIntent.begin(), addedIntent.end());
    *joins = rangeJoins(*joins, fuse);
  }
  return history.finish(store(kernel, fuse.Shape()));
}

EvolutionData filletRange(OcctKernel&                kernel,
                          uint32_t                   shapeId,
                          uint32_t                   edgeId,
                          double                     radius,
                          double                     start,
                          double                     end,
                          bool                       propagate,
                          const std::vector<int>&    hashes,
                          int                        bound,
                          const std::vector<double>& law)
{
  const char* op       = "filletRange";
  const auto& shape    = get(kernel, shapeId);
  const auto  selected = subShapes(kernel, op, shape, {edgeId}, TopAbs_EDGE);
  const auto  edge     = TopoDS::Edge(selected[0]);
  if (!std::isfinite(radius) || radius <= 0 || !std::isfinite(start) || !std::isfinite(end)
      || start < 0 || end > 1 || start >= end)
    fail(op, invalid_input, "positive radius and ordered range within [0,1] required");
  BRepAdaptor_Curve rangeCurve(edge);
  const double      rangeLength = GCPnts_AbscissaPoint::Length(
    rangeCurve,
    rangeCurve.FirstParameter()
      + (rangeCurve.LastParameter() - rangeCurve.FirstParameter()) * start,
    rangeCurve.FirstParameter() + (rangeCurve.LastParameter() - rangeCurve.FirstParameter()) * end);
  if (rangeLength <= BRepExt::tolerance::length)
    fail(op, invalid_input, "empty range: limits are within kernel resolution");
  BRepFilletAPI_MakeFillet maker(shape);
  maker.Add(radius, edge);
  const int contour = maker.Contour(edge);
  if (!law.empty())
  {
    if (law.size() < 4 || law.size() % 2 || law.front() != 0 || law[law.size() - 2] != 1)
      fail(op, invalid_input, "radius points must run from 0 to 1");
    const int size  = int(law.size() / 2);
    int       index = 1;
    while (index <= maker.NbEdges(contour) && !maker.Edge(contour, index).IsSame(edge))
      ++index;
    if (index > maker.NbEdges(contour))
      fail(op, no_solution, "selected edge is absent from contour");
    const bool reverse = maker.Edge(contour, index).Orientation() != edge.Orientation();
    NCollection_Array1<gp_Pnt2d> points(1, size);
    for (int i = 0; i < size; ++i)
    {
      if (!std::isfinite(law[2 * i]) || !std::isfinite(law[2 * i + 1]) || law[2 * i + 1] <= 0
          || (i && law[2 * i] <= law[2 * i - 2]))
        fail(op, invalid_input, "ordered finite parameters and positive radii required");
      points(reverse ? size - i : i + 1) =
        gp_Pnt2d(reverse ? 1 - law[2 * i] : law[2 * i], law[2 * i + 1]);
    }
    maker.SetRadius(points, contour, index);
  }
  const bool limited = start > 0 || end < 1 || (!propagate && maker.NbEdges(contour) > 1);
  maker.Build();
  if (!maker.IsDone() || !BRepExt::invalidity(maker.Shape(), shape).empty())
    fail(op, no_solution, "fillet construction failed");
  FaceHistory history(op, hashes, bound);
  if (!limited)
  {
    history.collect(maker, shape);
    return history.finish(store(kernel, maker.Shape()));
  }
  std::vector<TopoDS_Edge> others;
  for (int e = 1; e <= maker.NbEdges(contour); ++e)
    if (!maker.Edge(contour, e).IsSame(edge))
      others.push_back(maker.Edge(contour, e));
  return withinRange(kernel, op, shape, maker.Shape(), edge, start, end, others, history);
}

EvolutionData blendRange(OcctKernel&                  kernel,
                         uint32_t                     shapeId,
                         uint32_t                     edgeId,
                         double                       first,
                         double                       second,
                         blending::Law                law,
                         double                       rho,
                         uint32_t                     reference,
                         double                       angle,
                         double                       start,
                         double                       end,
                         const std::vector<int>&      hashes,
                         int                          bound,
                         const std::vector<double>&   sizes,
                         std::vector<blending::Join>* joins)
{
  const char* op    = "filletRange";
  const auto& shape = get(kernel, shapeId);
  const auto  edge  = TopoDS::Edge(subShapes(kernel, op, shape, {edgeId}, TopAbs_EDGE)[0]);
  if (!std::isfinite(start) || !std::isfinite(end) || start < 0 || end > 1 || start >= end)
    fail(op, invalid_input, "an ordered range within [0,1] required");
  if (law < blending::Law::Conic || law > blending::Law::Chamfer)
    fail(op, invalid_input, "a range takes a conic, G2 or chamfer section");
  if (sizes.size() % 2)
    fail(op, invalid_input, "size places come in (place, size) pairs");
  // The edge alone, no tangent edges added: the range lies on it.
  blending::Request request;
  request.operation = op;
  request.body      = shape;
  request.edges     = {edge};
  request.mayAdd    = {edge};
  request.first     = first;
  request.second    = second;
  request.law       = law;
  request.rho       = rho;
  request.angle     = angle;
  if (reference)
    request.reference = subShapes(kernel, op, shape, {reference}, TopAbs_FACE)[0];
  if (!sizes.empty())
  {
    auto& points = request.radii.emplace_back();
    for (size_t k = 0; k + 1 < sizes.size(); k += 2)
      points.emplace_back(sizes[k], sizes[k + 1]);
  }
  blending::Blended blended;
  try
  {
    blended = blending::blend(request);
  }
  catch (const Standard_Failure& error)
  {
    fail(op,
         no_solution,
         std::string("OpenCascade could not build the blend (") + error.ExceptionType() + ")");
  }
  double size = std::max(first, second);
  for (const auto& points : request.radii)
    for (const auto& [place, at] : points)
      size = std::max(size, std::max(first, second) * at / first);
  blending::accepted(op, shape, blended.shape, size, &blended);
  if (joins)
    *joins = blended.joins;
  FaceHistory history(op, hashes, bound);
  if (start == 0 && end == 1)
  {
    history.collect(*blended.history, shape);
    return history.finish(store(kernel, blended.shape));
  }
  return withinRange(kernel, op, shape, blended.shape, edge, start, end, {}, history, joins, &blended);
}
} // namespace BRepExt
