// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Shared by the TKAlgoExt tests: the host seam of BRepToolsExt_Session.hxx (a vector of shapes,
// defined in BRepExtTest_Support.cxx) and oracles that use nothing of the code under test.

#pragma once

#include <BRepToolsExt_Session.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pnt2d.hxx>

#include <cmath>
#include <functional>
#include <vector>

namespace BRepExtTest
{

OcctKernel& session();

/// The volume by adaptive Gauss-Kronrod integration to 1e-8 relative accuracy. The default fixed
/// rule is off by about 1e-5 of the volume on spline blend faces.
inline double volume(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::VolumePropertiesGK(shape, props, 1e-8, false, true);
  return props.Mass();
}

inline double length(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::LinearProperties(shape, props);
  return props.Mass();
}

inline double area(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::SurfaceProperties(shape, props);
  return props.Mass();
}

/// The sub-shapes of `shape` of one type, each once, whose bounding box `pick` accepts (the box
/// limits as xmin ymin zmin xmax ymax zmax). An explorer meets an edge once for every face and
/// wire that holds it.
inline std::vector<TopoDS_Shape> subShapesWhere(const TopoDS_Shape&                       shape,
                                                TopAbs_ShapeEnum                          type,
                                                const std::function<bool(const double*)>& pick)
{
  std::vector<TopoDS_Shape> found;
  for (TopExp_Explorer it(shape, type); it.More(); it.Next())
  {
    bool seen = false;
    for (const TopoDS_Shape& earlier : found)
      seen = seen || earlier.IsSame(it.Current());
    if (seen)
      continue;
    Bnd_Box box;
    BRepBndLib::Add(it.Current(), box);
    double b[6];
    box.Get(b[0], b[1], b[2], b[3], b[4], b[5]);
    if (pick(b))
      found.push_back(it.Current());
  }
  return found;
}

inline std::vector<TopoDS_Edge> edgesWhere(const TopoDS_Shape&                       shape,
                                           const std::function<bool(const double*)>& pick)
{
  std::vector<TopoDS_Edge> edges;
  for (const TopoDS_Shape& edge : subShapesWhere(shape, TopAbs_EDGE, pick))
    edges.push_back(TopoDS::Edge(edge));
  return edges;
}

inline std::vector<TopoDS_Face> facesWhere(const TopoDS_Shape&                       shape,
                                           const std::function<bool(const double*)>& pick)
{
  std::vector<TopoDS_Face> faces;
  for (const TopoDS_Shape& face : subShapesWhere(shape, TopAbs_FACE, pick))
    faces.push_back(TopoDS::Face(face));
  return faces;
}

inline bool near(double a, double b)
{
  return std::abs(a - b) < 1e-6;
}

struct Defect
{
  double worst = 0, rms = 0;
  int    samples = 0;
};

/// How far the blend faces of `result` are from the exact rolling-ball blend of radius r on
/// `input`, measured without help from the blend. A point P of the exact blend with unit normal N
/// has the ball centre C = P -/+ rN at distance exactly r from the input's faces, so the defect
/// dist(C, input) - r is zero for the exact blend, positive where the surface lies inside it (more
/// material removed) and negative outside. Faces lying on the input are not blends and are skipped;
/// each blend face takes the side of the normal with the smaller mean square. `rolls` is false
/// where the ball cannot roll (a centre inside a face the blend runs out onto). Distances are taken
/// to the faces of the input, because a solid is a volume to BRepExtrema_DistShapeShape and every
/// centre inside it would be at 0. 24 x 24 interior samples per face.
inline Defect ballDefect(
  const TopoDS_Shape&                       result,
  const TopoDS_Shape&                       input,
  double                                    r,
  const std::function<bool(const gp_Pnt&)>& rolls = [](const gp_Pnt&) { return true; })
{
  Defect          out;
  double          sum = 0;
  constexpr int   N   = 24;
  TopoDS_Compound skin;
  BRep_Builder    builder;
  builder.MakeCompound(skin);
  for (TopExp_Explorer f(input, TopAbs_FACE); f.More(); f.Next())
    builder.Add(skin, f.Current());
  for (TopExp_Explorer f(result, TopAbs_FACE); f.More(); f.Next())
  {
    const TopoDS_Face face = TopoDS::Face(f.Current());
    double            u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    const BRepAdaptor_Surface              surface(face);
    std::vector<std::pair<gp_Pnt, gp_Dir>> samples;
    for (int i = 1; i < N; ++i)
      for (int j = 1; j < N; ++j)
      {
        const double u = u0 + (u1 - u0) * i / N, v = v0 + (v1 - v0) * j / N;
        if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() != TopAbs_IN)
          continue;
        BRepLProp_SLProps props(surface, u, v, 1, 1e-9);
        if (!props.IsNormalDefined())
          continue;
        gp_Dir normal = props.Normal();
        if (face.Orientation() == TopAbs_REVERSED)
          normal.Reverse();
        samples.emplace_back(surface.Value(u, v), normal);
      }
    bool onInput = true;
    for (size_t k = 0; k < samples.size() && onInput; k += 7)
    {
      const BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(samples[k].first).Vertex(),
                                         skin);
      onInput = d.IsDone() && d.Value() <= 1e-6;
    }
    if (samples.empty() || onInput)
      continue;
    double bestSum = 0, bestWorst = 0;
    int    bestCount = 0;
    for (const double side : {-1., 1.})
    {
      double faceSum = 0, faceWorst = 0;
      int    count = 0;
      for (const auto& [point, normal] : samples)
      {
        const gp_Pnt centre = point.Translated(gp_Vec(normal) * (side * r));
        if (!rolls(centre))
          continue;
        const BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(centre).Vertex(), skin);
        if (!d.IsDone())
          continue;
        const double defect = d.Value() - r;
        faceSum += defect * defect;
        ++count;
        if (std::abs(defect) > std::abs(faceWorst))
          faceWorst = defect;
      }
      if (count && (bestCount == 0 || faceSum / count < bestSum / bestCount))
      {
        bestSum   = faceSum;
        bestWorst = faceWorst;
        bestCount = count;
      }
    }
    sum += bestSum;
    out.samples += bestCount;
    if (std::abs(bestWorst) > std::abs(out.worst))
      out.worst = bestWorst;
  }
  if (out.samples)
    out.rms = std::sqrt(sum / out.samples);
  return out;
}

/// The volume of a box a x b x c with every edge rounded by r (Steiner's formula for the parallel
/// body of the inner box): corners are eighths of a sphere, edges quarter cylinders, faces slabs.
inline double roundedBoxVolume(double a, double b, double c, double r)
{
  const double x = a - 2 * r, y = b - 2 * r, z = c - 2 * r;
  return x * y * z + 2 * r * (x * y + y * z + z * x) + M_PI * r * r * (x + y + z)
         + 4. / 3 * M_PI * r * r * r;
}

} // namespace BRepExtTest
