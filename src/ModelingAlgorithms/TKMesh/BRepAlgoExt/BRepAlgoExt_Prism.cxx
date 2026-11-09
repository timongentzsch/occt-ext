// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Resolution.hxx>
#include <BRepAlgoExt_Prism.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Trsf.hxx>
#include <cmath>
#include <algorithm>

namespace BRepExt
{
bool blendConsumesPrismCap(const TopoDS_Shape&              shape,
                           const std::vector<TopoDS_Shape>& edges,
                           double                           radius)
{
  // Limit this diagnostic to spline prisms. Polyhedral chamfers deliberately
  // permit a face to disappear through the existing supporting-plane fallback.
  bool spline = false;
  for (const auto& edge : edges)
    spline |= BRepAdaptor_Curve(TopoDS::Edge(edge)).GetType() == GeomAbs_BSplineCurve;
  if (!spline)
    return false;
  for (TopExp_Explorer faces(shape, TopAbs_FACE); faces.More(); faces.Next())
  {
    const auto          face = TopoDS::Face(faces.Current());
    BRepAdaptor_Surface cap(face);
    if (cap.GetType() != GeomAbs_Plane)
      continue;
    bool adjacent = false;
    for (TopExp_Explorer boundary(face, TopAbs_EDGE); boundary.More(); boundary.Next())
      for (const auto& edge : edges)
        adjacent |= edge.IsSame(boundary.Current());
    if (!adjacent)
      continue;
    const auto normal = cap.Plane().Axis().Direction();
    bool       prism  = true;
    int        caps   = 0;
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    {
      BRepAdaptor_Surface wall(TopoDS::Face(it.Current()));
      switch (wall.GetType())
      {
        case GeomAbs_Plane: {
          const double dot = std::abs(normal.Dot(wall.Plane().Axis().Direction()));
          if (dot > 1 - BRepExt::tolerance::fineAngle)
            ++caps;
          else if (dot > BRepExt::tolerance::fineAngle)
            prism = false;
          break;
        }
        case GeomAbs_SurfaceOfExtrusion:
          prism &= normal.IsParallel(wall.Direction(), BRepExt::tolerance::fineAngle);
          break;
        case GeomAbs_Cylinder:
          prism &=
            normal.IsParallel(wall.Cylinder().Axis().Direction(), BRepExt::tolerance::fineAngle);
          break;
        default:
          prism = false;
      }
    }
    if (!prism || caps != 2)
      continue;
    gp_Trsf toPlane;
    toPlane.SetTransformation(cap.Plane().Position());
    const auto flat = BRepBuilderAPI_Transform(face, toPlane, true).Shape();
    Bnd_Box    box;
    BRepBndLib::AddOptimal(flat, box, false, false);
    if (box.IsVoid() || box.IsOpen())
      continue;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    // Every disk contained by the cap must fit inside this enclosing strip.
    // Beyond half its width no inward-offset island can remain.
    if (2 * radius > std::min(xmax - xmin, ymax - ymin) + 2e-7)
      return true;
  }
  return false;
}
} // namespace BRepExt

#include <BRepTools_WireExplorer.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRep_Tool.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedDataMap.hxx>

namespace BRepExt
{
// Resolve the direction into a trimmed support in UV. This is independent of
// topological orientation, including a reflection and spline-converted planes.
bool supportInward(const TopoDS_Face& face,
                   const TopoDS_Edge& edge,
                   double             t,
                   const gp_Vec&      tangent,
                   gp_Vec&            inward)
{
  BRepAdaptor_Curve2d pc(edge, face);
  const auto          uv = pc.Value(t);
  BRepAdaptor_Surface surface(face);
  gp_Pnt              p;
  gp_Vec              du, dv;
  surface.D1(uv.X(), uv.Y(), p, du, dv);
  auto across = du.Crossed(dv).Crossed(tangent);
  if (across.SquareMagnitude() < 1e-20)
    return false;
  across.Normalize();
  const double aa = du.SquareMagnitude(), bb = dv.SquareMagnitude(), ab = du.Dot(dv),
               det = aa * bb - ab * ab;
  if (det < 1e-24)
    return false;
  const double u = (across.Dot(du) * bb - across.Dot(dv) * ab) / det;
  const double v = (across.Dot(dv) * aa - across.Dot(du) * ab) / det;
  for (double step : {1e-4, 1e-3, 1e-2})
    for (double sign : {1.0, -1.0})
    {
      BRepClass_FaceClassifier test(face,
                                    gp_Pnt2d(uv.X() + sign * step * u, uv.Y() + sign * step * v),
                                    1e-8);
      if (test.State() == TopAbs_IN)
      {
        inward = across * sign;
        return true;
      }
    }
  return false;
}

bool blendExceedsConvexRadius(const TopoDS_Shape&              shape,
                              const std::vector<TopoDS_Shape>& edges,
                              double                           radius,
                              const std::vector<TopoDS_Shape>* supports)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    adj;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adj);
  for (const auto& selected : edges)
  {
    const auto        edge = TopoDS::Edge(selected);
    BRepAdaptor_Curve curve(edge);
    if (!adj.Contains(edge) || adj.FindFromKey(edge).Extent() != 2)
      continue;
    const auto&  beside = adj.FindFromKey(edge);
    const double start = curve.FirstParameter(), span = curve.LastParameter() - start;
    const double middle = start + span * .5;
    gp_Pnt       point;
    gp_Vec       tangent;
    curve.D1(middle, point, tangent);
    for (const auto& item : beside)
    {
      const auto face = TopoDS::Face(item);
      if (supports && std::none_of(supports->begin(), supports->end(), [&](const auto& s) {
            return s.IsSame(face);
          }))
        continue;
      BRepAdaptor_Surface surface(face);
      if (surface.GetType() == GeomAbs_Plane)
        continue;
      const auto other = TopoDS::Face(beside.First().IsSame(face) ? beside.Last() : beside.First());
      gp_Vec     intoOther;
      if (!supportInward(other, edge, middle, tangent, intoOther))
        continue;
      BRepAdaptor_Curve2d pcurve(edge, face);
      const auto          mid = pcurve.Value(middle);
      gp_Pnt              p;
      gp_Vec              du, dv;
      surface.D1(mid.X(), mid.Y(), p, du, dv);
      const double side = intoOther.Dot(du.Crossed(dv));
      if (std::abs(side) < 1e-20)
        continue;
      // The ball lies on the side of this support toward the interior of
      // its neighbour, for convex and concave edges alike. The principal
      // bends are eigenvalues of II relative to I, not edge curvature or
      // cap setback. Evaluate the surface's exact derivatives.
      const int count =
        (surface.GetType() == GeomAbs_Cylinder || surface.GetType() == GeomAbs_Sphere) ? 1 : 4096;
      for (int i = 0; i < count; ++i)
      {
        const auto uv = pcurve.Value(start + span * (i + .5) / count);
        gp_Vec     uu, vv, uvd;
        surface.D2(uv.X(), uv.Y(), p, du, dv, uu, vv, uvd);
        auto normal = du.Crossed(dv);
        if (normal.SquareMagnitude() < 1e-24)
          continue;
        normal.Normalize();
        if (side < 0)
          normal.Reverse();
        const double E = du.SquareMagnitude(), F = du.Dot(dv), G = dv.SquareMagnitude();
        const double L = uu.Dot(normal), M = uvd.Dot(normal), N = vv.Dot(normal);
        const double a = E * G - F * F, b = E * N + G * L - 2 * F * M, c = L * N - M * M;
        if (a <= 1e-24)
          continue;
        const double bend = (b + std::sqrt(std::max(0.0, b * b - 4 * a * c))) / (2 * a);
        if (radius * bend > 1 + 1e-7)
          return true;
      }
    }
  }
  return false;
}
} // namespace BRepExt
