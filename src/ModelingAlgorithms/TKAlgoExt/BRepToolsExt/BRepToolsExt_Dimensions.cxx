// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>

namespace BRepExt
{
namespace
{
struct Primitive
{
  int    kind = 0;
  gp_Pnt origin;
  gp_Dir x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1};
  double a = 0, b = 0, c = 0;
};

Primitive recognize(const TopoDS_Shape& shape, double tolerance)
{
  Primitive                                                     p;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> solids, faces, vertices;
  TopExp::MapShapes(shape, TopAbs_SOLID, solids);
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
  if (solids.Extent() != 1)
    return p;
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  const double volume  = std::abs(props.Mass());
  auto         matches = [&](double expected) {
    return std::abs(volume - expected)
           <= std::max(tolerance * tolerance * tolerance, expected * 1e-8);
  };
  if (faces.Extent() == 1)
  {
    BRepAdaptor_Surface surface(TopoDS::Face(faces(1)));
    if (surface.GetType() == GeomAbs_Sphere)
    {
      const auto   sphere = surface.Sphere();
      const double r      = sphere.Radius();
      if (matches(4 * M_PI * r * r * r / 3))
      {
        p.kind   = 3;
        p.origin = sphere.Location();
        p.a      = 2 * r;
        return p;
      }
    }
  }
  if (faces.Extent() == 3)
  {
    int                 cylinders = 0, planes = 0;
    gp_Cylinder         cylinder;
    std::vector<gp_Pln> caps;
    for (int i = 1; i <= faces.Extent(); ++i)
    {
      BRepAdaptor_Surface s(TopoDS::Face(faces(i)));
      if (s.GetType() == GeomAbs_Cylinder)
      {
        cylinder = s.Cylinder();
        ++cylinders;
      }
      else if (s.GetType() == GeomAbs_Plane)
      {
        caps.push_back(s.Plane());
        ++planes;
      }
    }
    if (cylinders == 1 && planes == 2)
    {
      const auto axis = cylinder.Axis();
      const auto dir  = axis.Direction();
      if (std::abs(caps[0].Axis().Direction().Dot(dir)) < 1 - 1e-10
          || std::abs(caps[1].Axis().Direction().Dot(dir)) < 1 - 1e-10)
        return p;
      double a = gp_Vec(axis.Location(), caps[0].Location()).Dot(gp_Vec(dir)),
             b = gp_Vec(axis.Location(), caps[1].Location()).Dot(gp_Vec(dir));
      if (a > b)
        std::swap(a, b);
      if (matches(M_PI * cylinder.Radius() * cylinder.Radius() * (b - a)))
      {
        p.kind   = 2;
        p.origin = axis.Location().Translated(gp_Vec(dir) * a);
        p.z      = dir;
        p.x      = cylinder.XAxis().Direction();
        p.y      = cylinder.YAxis().Direction();
        p.a      = 2 * cylinder.Radius();
        p.b      = b - a;
        return p;
      }
    }
  }
  if (faces.Extent() != 6 || vertices.Extent() != 8)
    return p;
  std::vector<gp_Dir> axes;
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    BRepAdaptor_Surface s(TopoDS::Face(faces(i)));
    if (s.GetType() != GeomAbs_Plane)
      return p;
    auto n     = s.Plane().Axis().Direction();
    bool known = false;
    for (auto d : axes)
      if (std::abs(n.Dot(d)) > 1 - 1e-10)
        known = true;
    if (!known)
      axes.push_back(n);
  }
  if (axes.size() != 3)
    return p;
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
      if (std::abs(axes[i].Dot(axes[j])) > 1e-10)
        return p;
  // Deterministic world alignment, including boxes rotated away from world
  // axes.
  auto dominant = [](gp_Dir d) {
    return std::abs(d.X()) >= std::max(std::abs(d.Y()), std::abs(d.Z())) ? 0
           : std::abs(d.Y()) >= std::abs(d.Z())                          ? 1
                                                                         : 2;
  };
  std::sort(axes.begin(), axes.end(), [&](gp_Dir a, gp_Dir b) {
    return dominant(a) < dominant(b);
  });
  for (auto& d : axes)
    if (d.Coord(dominant(d) + 1) < 0)
      d.Reverse();
  axes[2]             = gp_Dir(gp_Vec(axes[0]).Crossed(gp_Vec(axes[1])));
  const gp_Pnt anchor = BRep_Tool::Pnt(TopoDS::Vertex(vertices(1)));
  double       lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
  for (int i = 1; i <= 8; ++i)
    for (int j = 0; j < 3; ++j)
    {
      double v = gp_Vec(anchor, BRep_Tool::Pnt(TopoDS::Vertex(vertices(i)))).Dot(gp_Vec(axes[j]));
      lo[j]    = std::min(lo[j], v);
      hi[j]    = std::max(hi[j], v);
    }
  for (int i = 1; i <= 8; ++i)
    for (int j = 0; j < 3; ++j)
    {
      double v = gp_Vec(anchor, BRep_Tool::Pnt(TopoDS::Vertex(vertices(i)))).Dot(gp_Vec(axes[j]));
      if (std::min(std::abs(v - lo[j]), std::abs(v - hi[j])) > tolerance)
        return p;
    }
  if (!matches((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])))
    return p;
  p.kind = 1;
  p.origin =
    anchor.Translated(gp_Vec(axes[0]) * lo[0] + gp_Vec(axes[1]) * lo[1] + gp_Vec(axes[2]) * lo[2]);
  p.x = axes[0];
  p.y = axes[1];
  p.z = axes[2];
  p.a = hi[0] - lo[0];
  p.b = hi[1] - lo[1];
  p.c = hi[2] - lo[2];
  return p;
}
} // namespace

std::vector<double> primitiveDimensions(OcctKernel& kernel, uint32_t shape, double tolerance)
{
  const auto& input = get(kernel, shape);
  if (input.ShapeType() == TopAbs_FACE)
  {
    BRepAdaptor_Surface surface(TopoDS::Face(input));
    if (surface.GetType() != GeomAbs_Cylinder
        || std::abs(surface.LastUParameter() - surface.FirstUParameter() - 2 * M_PI)
             > BRepExt::tolerance::angle)
      return {};
    const auto c = surface.Cylinder();
    const auto origin =
      c.Location().Translated(gp_Vec(c.Axis().Direction()) * surface.FirstVParameter());
    const auto x = c.XAxis().Direction(), y = c.YAxis().Direction(), z = c.Axis().Direction();
    return {2,
            origin.X(),
            origin.Y(),
            origin.Z(),
            x.X(),
            x.Y(),
            x.Z(),
            y.X(),
            y.Y(),
            y.Z(),
            z.X(),
            z.Y(),
            z.Z(),
            2 * c.Radius(),
            surface.LastVParameter() - surface.FirstVParameter(),
            0};
  }
  const auto p = recognize(input, tolerance);
  if (p.kind == 0)
    return {};
  return {double(p.kind),
          p.origin.X(),
          p.origin.Y(),
          p.origin.Z(),
          p.x.X(),
          p.x.Y(),
          p.x.Z(),
          p.y.X(),
          p.y.Y(),
          p.y.Z(),
          p.z.X(),
          p.z.Y(),
          p.z.Z(),
          p.a,
          p.b,
          p.c};
}

EvolutionData setPrimitiveDimensions(OcctKernel&             kernel,
                                     uint32_t                shape,
                                     int                     kind,
                                     double                  a,
                                     double                  b,
                                     double                  c,
                                     double                  tolerance,
                                     const std::vector<int>& hashes,
                                     int                     bound)
{
  const auto p = recognize(get(kernel, shape), tolerance);
  if (p.kind != kind || kind < 1 || kind > 3)
    fail("setPrimitiveDimensions",
         invalid_input,
         "body is not the requested unmodified primitive geometry");
  if (!std::isfinite(a) || a <= tolerance || !std::isfinite(b) || !std::isfinite(c)
      || (kind == 1 && (b <= tolerance || c <= tolerance)))
    fail("setPrimitiveDimensions", invalid_input, "dimensions must be positive");
  TopoDS_Shape result;
  gp_Ax2       axes(p.origin, p.z, p.x);
  if (kind == 1)
    result = BRepPrimAPI_MakeBox(axes, a, b, c).Shape();
  else if (kind == 2)
    result = BRepPrimAPI_MakeCylinder(axes, a / 2, p.b).Shape();
  else
    result = BRepPrimAPI_MakeSphere(p.origin, a / 2).Shape();
  if (kind == 1)
  {
    // A rectangular box has exactly one face for each outward axis direction.
    // Map those supports explicitly; keep analytic planes for subsequent edits.
    FaceHistory history("setPrimitiveDimensions", hashes, bound);
    auto        normal = [](const TopoDS_Shape& f) { return orientedPlaneAxis(TopoDS::Face(f)); };
    for (TopExp_Explorer old(get(kernel, shape), TopAbs_FACE); old.More(); old.Next())
    {
      if (!history.tracks(old.Current()))
        continue;
      NCollection_List<TopoDS_Shape> successors;
      for (TopExp_Explorer next(result, TopAbs_FACE); next.More(); next.Next())
        if (normal(old.Current()).Dot(normal(next.Current())) > 1 - 1e-10)
          successors.Append(next.Current());
      if (successors.Extent() != 1)
        fail("setPrimitiveDimensions", no_solution, "ambiguous box support evolution");
      history.modified(old.Current(), successors);
    }
    return history.finish(store(kernel, result));
  }
  return stored(kernel, result);
}
} // namespace BRepExt
