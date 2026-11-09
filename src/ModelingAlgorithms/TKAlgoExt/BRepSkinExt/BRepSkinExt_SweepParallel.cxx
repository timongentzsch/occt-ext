// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepSurfaceExt_FaceSeparation.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Law_Linear.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <gp_Pln.hxx>

namespace BRepExt
{
EvolutionData sweepSection(OcctKernel& kernel,
                           uint32_t    profileId,
                           uint32_t    pathId,
                           double      scale,
                           int         corner,
                           bool        parallel,
                           bool        solid,
                           double      tolerance)
{
  const char* op      = parallel ? "sweepParallel" : "sweepSection";
  const auto& profile = get(kernel, profileId);
  const auto& path    = get(kernel, pathId);
  if (profile.ShapeType() != TopAbs_WIRE || path.ShapeType() != TopAbs_WIRE || !std::isfinite(scale)
      || scale <= 0 || corner < 0 || corner > 2)
    fail(op,
         invalid_input,
         "planar profile and path wires, positive scale and valid corner mode "
         "required");
  if (!std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "tolerance must be positive");
  BRepOffsetAPI_MakePipeShell maker(TopoDS::Wire(path));
  if (parallel)
  {
    BRepBuilderAPI_MakeFace face(TopoDS::Wire(BRepBuilderAPI_Copy(profile).Shape()), true);
    gp_Ax2                  frame;
    if (face.IsDone())
    {
      BRepAdaptor_Surface surface(face.Face());
      if (surface.GetType() != GeomAbs_Plane)
        fail(op, invalid_input, "profile must be planar");
      frame = surface.Plane().Position().Ax2();
    }
    else
    {
      if (solid)
        fail(op, invalid_input, "solid profile must be planar and closed");
      // A single open straight edge has no unique supporting plane. Its fixed
      // section still translates along the spine, using the initial tangent.
      BRepAdaptor_CompCurve spine(TopoDS::Wire(path), true);
      gp_Pnt                at;
      gp_Vec                tangent;
      spine.D1(spine.FirstParameter(), at, tangent);
      if (tangent.Magnitude() <= BRepExt::tolerance::exactLength)
        fail(op, no_solution, "path has no initial tangent");
      frame = gp_Ax2(at, gp_Dir(tangent));
    }
    // A fixed section cannot reverse its signed flux through the spine. Include
    // every edge and its endpoints so a polyline corner cannot hide a reversal.
    int          sign = 0;
    const gp_Vec normal(frame.Direction());
    for (BRepTools_WireExplorer it(TopoDS::Wire(path)); it.More(); it.Next())
    {
      BRepAdaptor_Curve curve(it.Current());
      for (int i = 0; i <= 256; ++i)
      {
        gp_Pnt at;
        gp_Vec tangent;
        curve.D1(curve.FirstParameter()
                   + (curve.LastParameter() - curve.FirstParameter()) * i / 256,
                 at,
                 tangent);
        if (it.Current().Orientation() == TopAbs_REVERSED)
          tangent.Reverse();
        const double dot = tangent.Dot(normal);
        if (std::abs(dot) <= 1e-10 * tangent.Magnitude())
          fail(op, no_solution, "path tangent reaches the parallel section plane");
        const int current = dot > 0 ? 1 : -1;
        if (sign && current != sign)
          fail(op, no_solution, "path folds through the parallel section plane");
        sign = current;
      }
    }
    // Unlike SetMode(false) (corrected Frenet), this keeps every section parallel
    // to the authored profile plane, including on a curved spine.
    maker.SetMode(frame);
  }
  maker.SetTransitionMode(static_cast<BRepBuilderAPI_TransitionMode>(corner));
  maker.SetTolerance(tolerance, tolerance, BRepExt::tolerance::coarseAngle);
  if (scale == 1)
    maker.Add(profile, false, !parallel && corner != 1);
  else
  {
    BRepAdaptor_CompCurve curve(TopoDS::Wire(path), true);
    const double          length = GCPnts_AbscissaPoint::Length(curve);
    Handle(Law_Linear)    law    = new Law_Linear();
    law->Set(0, 1, length, scale);
    maker.SetLaw(profile, law, false, !parallel && corner != 1);
  }
  maker.Build();
  if (!maker.IsDone() || (solid && !maker.MakeSolid())
      || !BRepCheck_Analyzer(maker.Shape()).IsValid())
    fail(op, no_solution, "sweep did not produce the requested valid solid or sheet");
  // Only what the separation pre-pass cannot prove clear goes to the check.
  for (const auto& part : unprovenParts(maker.Shape()))
  {
    BOPAlgo_ArgumentAnalyzer analyzer;
    analyzer.SetShape1(part);
    analyzer.ArgumentTypeMode() = false;
    analyzer.SmallEdgeMode()    = false;
    analyzer.SelfInterMode()    = true;
    analyzer.Perform();
    if (analyzer.HasFaulty())
      fail(op, no_solution, "sweep intersects itself");
  }
  return stored(kernel, maker.Shape());
}
} // namespace BRepExt
