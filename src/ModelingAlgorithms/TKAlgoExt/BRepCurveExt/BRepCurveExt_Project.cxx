// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_PlanarCurve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>
#include <gp_Ax2.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>
#include <gp_Trsf.hxx>
#include <vector>

namespace BRepExt
{
namespace
{
gp_Dir direction(const char* op, Vec3 v)
{
  if (!finite(v) || gp_Vec(v.x, v.y, v.z).Magnitude() < 1e-15)
    fail(op, invalid_input, "direction must be finite and nonzero");
  return gp_Dir(v.x, v.y, v.z);
}

EvolutionData curves(OcctKernel& kernel, const char* op, const TopoDS_Shape& shape)
{
  if (shape.IsNull() || !TopExp_Explorer(shape, TopAbs_EDGE).More())
    fail(op, no_solution, "no curves were produced");
  if (!BRepCheck_Analyzer(shape).IsValid())
    fail(op, no_solution, "result curves are invalid");
  return stored(kernel, shape);
}

TopoDS_Shape curtain(OcctKernel& kernel, uint32_t wire, Vec3 dir, bool both)
{
  if (get(kernel, wire).ShapeType() != TopAbs_WIRE)
    fail("curveCurtain", invalid_input, "a curve wire is required");
  BRepPrimAPI_MakePrism prism(get(kernel, wire), direction("curveCurtain", dir), both, true);
  return prism.Shape();
}
} // namespace

EvolutionData curveCurtain(OcctKernel& kernel, uint32_t wire, Vec3 dir, bool both)
{
  return stored(kernel, curtain(kernel, wire, dir, both));
}

EvolutionData projectCurve(OcctKernel& kernel, uint32_t wire, uint32_t body, Vec3 dir, bool both)
{
  if (!TopExp_Explorer(get(kernel, body), TopAbs_FACE).More())
    fail("projectCurve", invalid_input, "target must contain faces");
  const char* op     = "projectCurve";
  const auto &source = get(kernel, wire), &target = get(kernel, body);
  if (source.ShapeType() != TopAbs_WIRE)
    fail(op, invalid_input, "source must be a wire");
  const gp_Vec travel(direction(op, dir));
  Bnd_Box      bounds;
  BRepBndLib::Add(source, bounds);
  BRepBndLib::Add(target, bounds);
  const double extent             = std::sqrt(bounds.SquareExtent()) * 2 + 1;
  TopoDS_Shape intersectionTarget = target;
  auto         sectionOf          = [&](TopoDS_Shape part) {
    if (both)
    {
      gp_Trsf shift;
      shift.SetTranslation(-travel * extent);
      part = BRepBuilderAPI_Transform(part, shift, true).Shape();
    }
    BRepPrimAPI_MakePrism prism(part, travel * (both ? 2 * extent : extent), true);
    BRepAlgoAPI_Section   section(prism.Shape(), intersectionTarget, false);
    section.Approximation(true);
    section.Build();
    if (!section.IsDone())
      fail(op, no_solution, "projection intersection failed");
    return section.Shape();
  };
  auto                           result = sectionOf(source);
  IntCurvesFace_ShapeIntersector rays;
  rays.Load(target, BRepExt::tolerance::length);
  auto complete = [&](const TopoDS_Shape& candidate) {
    DistanceTo toCandidate(candidate);
    for (TopExp_Explorer it(source, TopAbs_EDGE); it.More(); it.Next())
    {
      BRepAdaptor_Curve c(TopoDS::Edge(it.Current()));
      for (int i = 0; i <= 16; ++i)
      {
        const auto p =
          c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 16);
        rays.Perform(gp_Lin(p, gp_Dir(travel)), both ? -extent : 0, extent);
        for (int j = 1; j <= rays.NbPnt(); ++j)
        {
          const auto apart = toCandidate.from(rays.Pnt(j));
          if (!apart || *apart > 1e-5)
            return false;
        }
      }
    }
    return true;
  };
  if (!complete(result))
  {
    // A complete sphere has no geometric pole, only a singular UV chart.
    // Move that chart's axis away from the projection direction and retry the
    // exact same sphere locus, rather than approximating the missing branch.
    TopExp_Explorer it(target, TopAbs_FACE);
    auto            face = TopoDS::Face(it.Current());
    it.Next();
    BRepAdaptor_Surface surface(face);
    GProp_GProps        area;
    BRepGProp::SurfaceProperties(target, area);
    if (!it.More() && surface.GetType() == GeomAbs_Sphere
        && std::abs(area.Mass() - 4 * M_PI * std::pow(surface.Sphere().Radius(), 2)) < 1e-5)
    {
      const auto   sphere = surface.Sphere();
      const gp_Vec ref    = std::abs(travel.X()) < 0.8 ? gp_Vec(1, 0, 0) : gp_Vec(0, 1, 0);
      const gp_Dir axis(travel.Crossed(ref));
      intersectionTarget =
        BRepPrimAPI_MakeSphere(gp_Ax2(sphere.Location(), axis), sphere.Radius()).Shape();
      result = sectionOf(source);
    }
    if (!complete(result))
      fail(op, no_solution, "projection does not cover the source rays");
  }
  return curves(kernel, op, result);
}

EvolutionData projectPlane(OcctKernel& kernel, uint32_t shapeId, Vec3 origin, Vec3 normal)
{
  const auto         n     = direction("projectPlane", normal);
  Handle(Geom_Plane) plane = new Geom_Plane(gp_Pnt(origin.x, origin.y, origin.z), n);
  BRep_Builder       builder;
  TopoDS_Compound    result;
  builder.MakeCompound(result);
  std::vector<TopoDS_Edge> accepted;
  for (TopExp_Explorer it(get(kernel, shapeId), TopAbs_EDGE); it.More(); it.Next())
  {
    const auto edge = TopoDS::Edge(it.Current());
    double     a, b;
    const auto source = BRep_Tool::Curve(edge, a, b);
    if (source.IsNull())
      continue;
    const auto projected =
      GeomProjLib::ProjectOnPlane(new Geom_TrimmedCurve(source, a, b), plane, n, true);
    if (projected.IsNull())
      continue;
    BRepBuilderAPI_MakeEdge make(projected, a, b);
    if (!make.IsDone())
      continue; // An edge parallel to the projection direction collapses.
    auto         out = make.Edge();
    GProp_GProps props;
    BRepGProp::LinearProperties(out, props);
    if (props.Mass() <= BRepExt::tolerance::length)
      continue;
    // Coincident 3D silhouette edges can have different depths but the same
    // image.
    BRepAdaptor_Curve c(out);
    bool              duplicate = false;
    for (const auto& other : accepted)
    {
      Bnd_Box aBox, bBox;
      BRepBndLib::Add(out, aBox);
      BRepBndLib::Add(other, bBox);
      if (aBox.IsOut(bBox))
        continue;
      GProp_GProps q;
      BRepGProp::LinearProperties(other, q);
      if (std::abs(q.Mass() - props.Mass()) > BRepExt::tolerance::length)
        continue;
      bool same = true;
      for (int i = 0; i <= 4 && same; ++i)
      {
        auto p = c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 4);
        BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(p).Vertex(), other);
        same = d.IsDone() && d.Value() <= BRepExt::tolerance::length;
      }
      if (same)
      {
        duplicate = true;
        break;
      }
    }
    if (duplicate)
      continue;
    accepted.push_back(out);
    out.Orientation(edge.Orientation());
    builder.Add(result, out);
  }
  return curves(kernel, "projectPlane", result);
}

EvolutionData slot(OcctKernel& kernel, uint32_t wireId, double width, Vec3 normal)
{
  const auto& shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE || shape.Closed() || !std::isfinite(width) || width <= 0)
    fail("slot", invalid_input, "slot requires an open planar wire and positive width");
  // OCCT's closed-result offset of an open wire offsets both sides and creates
  // circular caps.
  TopExp_Explorer vertex(shape, TopAbs_VERTEX);
  if (!vertex.More())
    fail("slot", invalid_input, "curve has no vertices");
  const auto plane =
    gp_Pln(BRep_Tool::Pnt(TopoDS::Vertex(vertex.Current())), direction("slot", normal));
  for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next())
  {
    double     a, b;
    const auto curve = BRep_Tool::Curve(TopoDS::Edge(it.Current()), a, b);
    requirePlanarCurve("slot", curve, a, b, plane, BRepExt::tolerance::sewing);
  }
  BRepBuilderAPI_MakeFace  support(plane, TopoDS::Wire(shape), true);
  BRepOffsetAPI_MakeOffset offset(support.Face(), GeomAbs_Arc, false);
  offset.Perform(width / 2);
  if (!offset.IsDone())
    fail("slot", no_solution, "offset failed (nonplanar or self-intersecting input)");
  ShapeUpgrade_UnifySameDomain unify(offset.Shape(), true, false, false);
  unify.Build();
  return curves(kernel, "slot", unify.Shape());
}
} // namespace BRepExt
