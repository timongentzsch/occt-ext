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

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_HiddenLine.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBndLib.hxx>
#include <BRepLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pln.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <Geom_Plane.hxx>
#include <HLRBRep_Algo.hxx>
#include <HLRBRep_HLRToShape.hxx>
#include <HLRAlgo_Projector.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt
{
namespace
{
constexpr const char* op = "hiddenLine";

// All hatch families are curves on the exact trimmed faces. HLR then hides them with the
// same faces as the body's edges; a viewport mesh never decides visibility.
void hatchPlane(const TopoDS_Shape& shape,
                gp_Dir              axis,
                double              spacing,
                BRep_Builder&       builder,
                TopoDS_Compound&    hatch)
{
  if (spacing <= 0)
    return;
  Bnd_Box box;
  BRepBndLib::Add(shape, box, false);
  double a, b, c, d, e, f;
  box.Get(a, b, c, d, e, f);
  double lo = INFINITY, hi = -INFINITY;
  for (double x : {a, d})
    for (double y : {b, e})
      for (double z : {c, f})
      {
        const double at = gp_Vec(x, y, z).Dot(gp_Vec(axis));
        lo              = std::min(lo, at);
        hi              = std::max(hi, at);
      }
  if ((hi - lo) / spacing > 2000)
    fail(op, invalid_input, "hatching exceeds 2000 planes; increase spacing");
  for (double at = (std::floor(lo / spacing) + 1) * spacing; at < hi - 1e-6; at += spacing)
  {
    if (at <= lo + 1e-6)
      continue;
    BRepAlgoAPI_Section cut;
    cut.Init1(shape);
    cut.Init2(gp_Pln(gp_Pnt(axis.X() * at, axis.Y() * at, axis.Z() * at), axis));
    cut.SetNonDestructive(true);
    cut.Approximation(true);
    cut.Build();
    if (!cut.IsDone())
      fail(op, no_solution, "hatching plane intersection failed");
    for (TopExp_Explorer it(cut.Shape(), TopAbs_EDGE); it.More(); it.Next())
      builder.Add(hatch, it.Current());
  }
}

void isoHatch(const TopoDS_Face& face,
              bool               u,
              double             spacing,
              BRep_Builder&      builder,
              TopoDS_Compound&   hatch)
{
  if (spacing <= 0)
    return;
  if (1 / spacing > 2000)
    fail(op, invalid_input, "hatching exceeds 2000 isoparameters; increase spacing");
  const int count = int(std::ceil(1 / spacing));
  double    u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  const auto surface = BRep_Tool::Surface(face);
  for (int i = 1; i < count; ++i)
  {
    const double t     = i * spacing;
    const auto   curve = u ? surface->UIso(u0 + t * (u1 - u0)) : surface->VIso(v0 + t * (v1 - v0));
    BRepBuilderAPI_MakeEdge edge(curve, u ? v0 : u0, u ? v1 : u1);
    if (!edge.IsDone())
      continue;
    BRepAlgoAPI_Common clip;
    clip.SetArguments({edge.Shape()});
    clip.SetTools({face});
    clip.SetNonDestructive(true);
    clip.Build();
    if (!clip.IsDone())
      fail(op, no_solution, "isoparameter trim failed");
    for (TopExp_Explorer it(clip.Shape(), TopAbs_EDGE); it.More(); it.Next())
      builder.Add(hatch, it.Current());
  }
}
} // namespace

std::vector<DrawingCurve> hiddenLine(OcctKernel&                  kernel,
                                     const std::vector<uint32_t>& ids,
                                     Vec3                         normal,
                                     Vec3                         x,
                                     const std::vector<double>&   p)
{
  // chord, angle, smooth angle, planar spacing/axis, section spacing/axis,
  // radial around/along/about, U/V, section enabled/origin/normal.
  if ((p.size() != 24 && p.size() != 28) || ids.empty())
    fail(op, invalid_input, "drawing options and bodies are required");
  for (double value : p)
    if (!std::isfinite(value))
      fail(op, invalid_input, "drawing options must be finite");
  if (p[0] <= 0 || p[1] <= 0 || p[2] < 0)
    fail(op, invalid_input, "positive chord and angle tolerances are required");
  gp_Ax2 axes(p.size() == 28 ? gp_Pnt(p[24], p[25], p[26]) : gp_Pnt(0, 0, 0),
              gp_Dir(normal.x, normal.y, normal.z),
              gp_Dir(x.x, x.y, x.z));
  Handle(HLRBRep_Algo) algo = new HLRBRep_Algo;
  algo->Projector(p.size() == 28 && p[27] > 0 ? HLRAlgo_Projector(axes, p[27])
                                              : HLRAlgo_Projector(axes));
  std::vector<TopoDS_Shape> bodies, hatches, boundaries;
  BRep_Builder              builder;
  for (auto id : ids)
  {
    // HLR writes p-curves into its input. Work on a private graph, including on errors.
    TopoDS_Shape shape = BRepBuilderAPI_Copy(get(kernel, id), true, false).Shape();
    if (p[17] != 0)
    {
      const gp_Pnt    origin(p[18], p[19], p[20]);
      const gp_Dir    n(p[21], p[22], p[23]);
      const auto      plane = BRepBuilderAPI_MakeFace(gp_Pln(origin, n)).Face();
      const auto      half = BRepPrimAPI_MakeHalfSpace(plane, origin.Translated(gp_Vec(n))).Solid();
      BRepAlgoAPI_Cut cut;
      cut.SetArguments({shape});
      cut.SetTools({half});
      cut.SetNonDestructive(true);
      cut.Build();
      if (!cut.IsDone())
        fail(op, no_solution, "section cut failed");
      shape = cut.Shape();
    }
    if (p[2] > 0)
      BRepLib::EncodeRegularity(shape, p[2]);
    bodies.push_back(shape);
    TopoDS_Shape boundary;
    if (p[17] && TopExp_Explorer(shape, TopAbs_FACE).More())
    {
      // Reuse the exact outline for the sectioned body. HLR hides it against every body.
      const auto temporary = store(kernel, shape);
      try
      {
        const auto result = outline(kernel, temporary, {0, 0, 0}, normal, x, true);
        boundary          = get(kernel, result.resultId);
        release(kernel, result.resultId);
        release(kernel, temporary);
      }
      catch (...)
      {
        release(kernel, temporary);
        throw;
      }
    }
    boundaries.push_back(boundary);
    TopoDS_Compound hatch;
    builder.MakeCompound(hatch);
    if (!TopExp_Explorer(shape, TopAbs_EDGE).More())
    {
      hatches.push_back(hatch);
      continue;
    }
    if (p[3] > 0)
      hatchPlane(shape, gp_Dir(p[4], p[5], p[6]), p[3], builder, hatch);
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    {
      const auto                face = TopoDS::Face(it.Current());
      const BRepAdaptor_Surface surface(face);
      isoHatch(face, true, p[15], builder, hatch);
      isoHatch(face, false, p[16], builder, hatch);
      const auto type = surface.GetType();
      if (type == GeomAbs_Cylinder || type == GeomAbs_Sphere)
      {
        isoHatch(face, true, p[12], builder, hatch);
        if (type == GeomAbs_Cylinder)
          isoHatch(face, false, p[13], builder, hatch);
        if (type == GeomAbs_Sphere)
          isoHatch(face, false, p[14], builder, hatch);
      }
      if (p[17] && p[7] > 0 && type == GeomAbs_Plane)
      {
        const auto   plane = surface.Plane();
        const gp_Pnt origin(p[18], p[19], p[20]);
        const gp_Dir n(p[21], p[22], p[23]);
        if (plane.Distance(origin) < BRepExt::tolerance::sewing
            && plane.Axis().Direction().IsParallel(n, 1e-8))
          hatchPlane(face, gp_Dir(p[8], p[9], p[10]), p[7], builder, hatch);
      }
    }
    hatches.push_back(hatch);
    algo->Add(shape);
    if (!boundary.IsNull())
      algo->Add(boundary);
    if (TopExp_Explorer(hatch, TopAbs_EDGE).More())
      algo->Add(hatch);
  }
  if (std::none_of(bodies.begin(), bodies.end(), [](const TopoDS_Shape& s) {
        return TopExp_Explorer(s, TopAbs_EDGE).More();
      }))
    return {};
  algo->Update();
  algo->Hide();
  HLRBRep_HLRToShape        read(algo);
  std::vector<DrawingCurve> out;
  auto                      append = [&](const TopoDS_Shape& shape, int body, int kind) {
    if (shape.IsNull())
      return;
    for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next())
    {
      const auto edge = TopoDS::Edge(it.Current());
      if (BRep_Tool::Degenerated(edge))
        continue;
      BRepAdaptor_Curve           curve(edge);
      GCPnts_TangentialDeflection poly(curve, p[1], p[0], 2);
      DrawingCurve                line{body, kind, {}};
      for (int i = 1; i <= poly.NbPoints(); ++i)
      {
        const auto q = poly.Value(i);
        line.points.push_back(q.X());
        line.points.push_back(q.Y());
      }
      if (line.points.size() >= 4)
        out.push_back(std::move(line));
    }
  };
  for (size_t i = 0; i < bodies.size(); ++i)
  {
    const auto& shape = bodies[i];
    if (!TopExp_Explorer(shape, TopAbs_EDGE).More())
      continue;
    append(read.VCompound(shape), int(i), 0);
    append(read.HCompound(shape), int(i), 1);
    append(read.Rg1LineVCompound(shape), int(i), 2);
    append(read.RgNLineVCompound(shape), int(i), 2);
    append(read.Rg1LineHCompound(shape), int(i), 3);
    append(read.RgNLineHCompound(shape), int(i), 3);
    append(read.OutLineVCompound(shape), int(i), 4);
    append(read.OutLineHCompound(shape), int(i), 5);
    if (!boundaries[i].IsNull())
    {
      append(read.VCompound(boundaries[i]), int(i), 8);
      append(read.HCompound(boundaries[i]), int(i), 5);
    }
    if (TopExp_Explorer(hatches[i], TopAbs_EDGE).More())
    {
      append(read.VCompound(hatches[i]), int(i), 6);
      append(read.HCompound(hatches[i]), int(i), 7);
    }
  }
  return out;
}
} // namespace BRepExt
