// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Draft a planar profile directly, without independently cutting/fusing its adjoining walls.
#include <BRepBooleanExt_Accept.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepOffsetExt_SewHistory.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepOffsetAPI_MakeDraft.hxx>
#include <BRepTools.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Plane.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Lin.hxx>
#include <gp_Pnt2d.hxx>

namespace BRepExt
{
TopoDS_Shape draftExtrude(const TopoDS_Shape&        profile,
                          const gp_Vec&              vector,
                          double                     angle,
                          Handle(BRepTools_History)& history)
{
  constexpr const char* op = "extrude";
  if (!std::isfinite(angle) || std::abs(angle) <= Precision::Angular()
      || std::abs(angle) >= M_PI / 2)
    fail(op, invalid_input, "Draft angle must be finite, nonzero and below 90 degrees");
  if (profile.ShapeType() != TopAbs_FACE)
    fail(op, unsupported_geometry, "Direct draft requires a planar closed face");
  int         wires = 0;
  TopoDS_Wire wire;
  for (TopExp_Explorer it(profile, TopAbs_WIRE); it.More(); it.Next())
  {
    ++wires;
    wire = TopoDS::Wire(it.Current());
    if (!BRep_Tool::IsClosed(it.Current()))
      fail(op, unsupported_geometry, "Direct draft requires a closed profile wire");
  }
  // OCCT's face constructor uses only its first wire; never fill an omitted hole.
  if (wires != 1)
    fail(op, unsupported_geometry, "Direct draft requires exactly one profile wire");
  const BRepAdaptor_Surface support(TopoDS::Face(profile));
  const gp_Dir              pull(vector);
  if (support.GetType() != GeomAbs_Plane
      || std::abs(support.Plane().Axis().Direction().Dot(pull)) < 1 - tolerance::angle)
    fail(op,
         unsupported_geometry,
         "Direct draft requires a planar profile perpendicular to extrusion");
  if (!BRepCheck_Analyzer(profile).IsValid())
    fail(op, invalid_input, "The draft profile is invalid");
  // General curved boundaries use the existing draft path: MakeDraft's approximate
  // offset sweep is not guaranteed to meet the wall-angle certificate for them.
  for (TopExp_Explorer it(profile, TopAbs_EDGE); it.More(); it.Next())
  {
    const BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
    if (curve.GetType() == GeomAbs_Line || curve.GetType() == GeomAbs_Circle)
      continue;
    std::vector<gp_Pnt> poles;
    if (curve.GetType() == GeomAbs_BezierCurve)
    {
      const auto bezier = curve.Bezier();
      for (int i = 1; i <= bezier->NbPoles(); ++i)
        poles.push_back(bezier->Pole(i));
    }
    else if (curve.GetType() == GeomAbs_BSplineCurve)
    {
      const auto spline = curve.BSpline();
      for (int i = 1; i <= spline->NbPoles(); ++i)
        poles.push_back(spline->Pole(i));
    }
    if (poles.size() < 2 || poles.front().Distance(poles.back()) <= Precision::Confusion())
      fail(op, unsupported_geometry, "Direct draft requires line or circle profile edges");
    const gp_Lin line(poles.front(), gp_Dir(gp_Vec(poles.front(), poles.back())));
    for (const auto& p : poles)
      if (line.Distance(p) > Precision::Confusion())
        fail(op, unsupported_geometry, "Direct draft requires line or circle profile edges");
  }
  const auto              base   = support.Plane().Location();
  const double            height = vector.Magnitude();
  BRepOffsetAPI_MakeDraft maker(wire, pull, std::abs(angle));
  const gp_Dir            normal = orientedPlaneAxis(TopoDS::Face(profile));
  // Angle is relative to the oriented source face, not the unsigned pull: inward
  // positive draft widens the subtractive tool, matching the outward wall inclination.
  maker.SetDraft((angle > 0) == (normal.Dot(pull) > 0));
  // Perform(length) limits slanted generators, not the axial extrusion. Stop on the exact plane.
  const Handle(Geom_Surface) stop = new Geom_Plane(base.Translated(vector), pull);
  maker.Perform(stop, true);
  checkCancelled(op);
  if (!maker.IsDone())
    fail(op, no_solution, "The profile draft did not complete");
  // Keep the stationary cap out of MakeDraft's opaque internal sewing so its exact
  // edge/vertex correspondence remains available through our sewing history.
  BRepBuilderAPI_Sewing sewing(tolerance::sewing);
  sewing.Add(profile);
  for (TopExp_Explorer face(maker.Shape(), TopAbs_FACE); face.More(); face.Next())
    sewing.Add(face.Current());
  sewing.Perform();
  const auto result = closedSolid(op, sewing, invalid_result, 0);
  if (result.IsNull())
    fail(op, unsupported_geometry, "The profile draft cannot close this boundary as a solid");
  if (!accepts(result, true))
    fail(op, invalid_result, "The profile draft is invalid");

  history                        = new BRepTools_History();
  int                      bases = 0, tops = 0;
  TopoDS_Face              bottom;
  std::vector<TopoDS_Face> sides;
  for (TopExp_Explorer it(result, TopAbs_FACE); it.More(); it.Next())
  {
    checkCancelled(op);
    const auto          face = TopoDS::Face(it.Current());
    BRepAdaptor_Surface surface(face);
    if (surface.GetType() == GeomAbs_Plane
        && std::abs(surface.Plane().Axis().Direction().Dot(pull)) > 1 - tolerance::angle)
    {
      const double z = gp_Vec(base, surface.Plane().Location()).Dot(gp_Vec(pull));
      if (std::abs(z) <= tolerance::sewing)
      {
        ++bases;
        bottom = face;
        history->AddModified(profile, face);
      }
      else if (std::abs(z - height) <= tolerance::sewing)
      {
        ++tops;
        history->AddGenerated(profile, face);
      }
      else
        fail(op, invalid_result, "The profile draft has an unexpected cap");
      continue;
    }
    sides.push_back(face);
    // Verify the actual trimmed faces, not only the sweep's intermediate surfaces.
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    int samples = 0;
    for (int i = 0; i < 9; i++)
      for (int j = 0; j < 5; j++)
      {
        const double u = u0 + (u1 - u0) * (i + .5) / 9, v = v0 + (v1 - v0) * (j + .5) / 5;
        BRepClass_FaceClassifier inside(face, gp_Pnt2d(u, v), tolerance::sewing);
        if (inside.State() != TopAbs_IN && inside.State() != TopAbs_ON)
          continue;
        gp_Pnt at;
        gp_Vec du, dv;
        surface.D1(u, v, at, du, dv);
        const auto n = du.Crossed(dv);
        if (!std::isfinite(n.Magnitude()) || n.Magnitude() <= gp::Resolution())
          fail(op, invalid_result, "The profile draft has a singular wall");
        const double measured =
          std::asin(std::clamp(std::abs(n.Dot(gp_Vec(pull))) / n.Magnitude(), 0., 1.));
        const double z = gp_Vec(base, at).Dot(gp_Vec(pull));
        if (std::abs(measured - std::abs(angle)) > tolerance::angle || z < -tolerance::sewing
            || z > height + tolerance::sewing)
          fail(op,
               invalid_result,
               "The profile draft does not meet its requested angle and extent");
        ++samples;
      }
    if (!samples)
      fail(op, invalid_result, "The profile draft wall could not be checked");
  }
  if (bases != 1 || tops != 1 || sides.empty())
    fail(op, invalid_result, "The profile draft did not preserve its two caps");
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> baseEdges;
  TopExp::MapShapes(bottom, TopAbs_EDGE, baseEdges);
  if (!sewingHistory(sewing, {profile}, result, {TopAbs_EDGE, TopAbs_VERTEX}, *history).empty())
    fail(op, invalid_result, "The stationary draft boundary has incomplete sewing history");
  // Wire-input MakeDraft reports the terminal-cut wall faces. Carry each through our
  // explicit sewing history, requiring complete unambiguous ancestry for every wall.
  std::vector<int> owners(sides.size());
  for (TopExp_Explorer edge(profile, TopAbs_EDGE); edge.More(); edge.Next())
  {
    for (const auto& generated : maker.Generated(edge.Current()))
    {
      if (generated.ShapeType() != TopAbs_FACE)
        continue;
      const auto successor = sewn(sewing, generated);
      for (size_t i = 0; i < sides.size(); ++i)
      {
        if (!successor.IsSame(sides[i]))
          continue;
        ++owners[i];
        history->AddGenerated(edge.Current(), sides[i]);
        for (TopExp_Explorer boundary(sides[i], TopAbs_EDGE); boundary.More(); boundary.Next())
          if (!baseEdges.Contains(boundary.Current()))
            history->AddGenerated(edge.Current(), boundary.Current());
      }
    }
  }
  for (const int count : owners)
    if (count != 1)
      fail(op, invalid_result, "The profile draft has ambiguous generated-wall history");
  return result;
}
} // namespace BRepExt
