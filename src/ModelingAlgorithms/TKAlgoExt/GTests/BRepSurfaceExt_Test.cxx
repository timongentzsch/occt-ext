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

// Non-uniform scaling of BRepSurfaceExt on bodies with curved faces: the volume is multiplied by
// the product of the factors and the bounding box follows the factors exactly. Bridging and
// alignment of surfaces against the definition of G1: the normals of the two sides agree along the
// shared boundary, and the proofs of face separation against distances.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <gp_Ax1.hxx>
#include <limits>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepSurfaceExt_FaceSeparation.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

// Scaling x, y, z by (fx, fy, fz) about a centre multiplies every volume by fx fy fz and every
// point of every face lands where the factors put it. A sphere becomes an ellipsoid of volume 4/3
// pi fx fy fz R^3, a cylinder an elliptic one: faces that are no longer analytic must still be
// valid.
TEST(BRepSurfaceExt, NonUniformScaleMultipliesVolumeByTheProductOfFactors)
{
  const double        R = 6;
  const BRepExt::Vec3 factors{2.0, 0.5, 3.0}, centre{1.0, -2.0, 4.0};
  const TopoDS_Shape  bodies[2] = {
    BRepPrimAPI_MakeSphere(gp_Pnt(1, -2, 4), R).Shape(),
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, 10).Shape()};
  for (const TopoDS_Shape& body : bodies)
  {
    const auto          result = BRepExt::scaleAxes(session(),
                                           BRepExt::store(session(), body),
                                           factors,
                                           centre,
                                                    {},
                                           BRepExt::failureHashBound);
    const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made),
                volume(body) * factors.x * factors.y * factors.z,
                1e-8 * volume(body) * 3);
    // Every point of the result, mapped back by the inverse of the scaling, lies on the original
    // body's skin: the faces were scaled, not approximated.
    TopoDS_Compound skin;
    BRep_Builder    builder;
    builder.MakeCompound(skin);
    for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
      builder.Add(skin, f.Current());
    double worst = 0;
    for (TopExp_Explorer f(made, TopAbs_FACE); f.More(); f.Next())
    {
      const TopoDS_Face face = TopoDS::Face(f.Current());
      double            u0, u1, v0, v1;
      BRepTools::UVBounds(face, u0, u1, v0, v1);
      const BRepAdaptor_Surface surface(face);
      for (int i = 1; i < 12; ++i)
        for (int j = 1; j < 12; ++j)
        {
          const double u = u0 + (u1 - u0) * i / 12, v = v0 + (v1 - v0) * j / 12;
          if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() != TopAbs_IN)
            continue;
          const gp_Pnt                     p = surface.Value(u, v);
          const gp_Pnt                     original(centre.x + (p.X() - centre.x) / factors.x,
                                centre.y + (p.Y() - centre.y) / factors.y,
                                centre.z + (p.Z() - centre.z) / factors.z);
          const BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(original).Vertex(), skin);
          ASSERT_TRUE(d.IsDone());
          worst = std::max(worst, d.Value());
        }
    }
    EXPECT_LE(worst, 1e-7);
  }
}

namespace
{

TopoDS_Face planeFace(const gp_Pnt& origin,
                      const gp_Dir& normal,
                      const gp_Dir& xDirection,
                      double        u0,
                      double        u1,
                      double        v0,
                      double        v1)
{
  return BRepBuilderAPI_MakeFace(new Geom_Plane(gp_Ax3(origin, normal, xDirection)),
                                 u0,
                                 u1,
                                 v0,
                                 v1,
                                 1e-7)
    .Face();
}

//! The unit normal of the face at its parameter (u, v), with the face's own orientation ignored.
gp_Vec normalAt(const TopoDS_Face& face, double u, double v)
{
  BRepAdaptor_Surface surface(face);
  gp_Pnt              p;
  gp_Vec              du, dv;
  surface.D1(u, v, p, du, dv);
  return du.Crossed(dv).Normalized();
}

} // namespace

// A G1 bridge between the edge of a floor and the edge of a wall, both along y: the floor z = 0
// ends at x = 10, the wall x = 16 begins at z = 3 and rises. The bridge must start on the floor's
// edge with the floor's normal (0, 0, 1) and arrive on the wall's edge with the wall's normal
// (1, 0, 0), up to sign, and stay inside the box between the two edges, which holds for a Hermite
// blend because its control polygon does. A bridge that is only G0 meets both edges but not these
// normals; the tangent of a wrong end shows as an angle of order one.
TEST(BRepSurfaceExt, G1BridgeBetweenFloorAndWallHasTheNormalsOfBothAtTheirEdges)
{
  const TopoDS_Face floor =
    planeFace(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0), 0, 10, 0, 6);
  const TopoDS_Face wall =
    planeFace(gp_Pnt(16, 0, 3), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0), 0, 6, 0, 10);
  const auto floorEdge =
    edgesWhere(floor, [](const double* l) { return near(l[0], 10) && near(l[3], 10); });
  const auto wallEdge =
    edgesWhere(wall, [](const double* l) { return near(l[2], 3) && near(l[5], 3); });
  ASSERT_EQ(floorEdge.size(), 1u);
  ASSERT_EQ(wallEdge.size(), 1u);
  const auto          result = BRepExt::bridgeSurface(session(),
                                             BRepExt::store(session(), floorEdge[0]),
                                             BRepExt::store(session(), floor),
                                             BRepExt::store(session(), wallEdge[0]),
                                             BRepExt::store(session(), wall),
                                             1,
                                             1,
                                             1.0,
                                             1.0,
                                             1e-7,
                                             0);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  TopExp_Explorer     it(made, TopAbs_FACE);
  ASSERT_TRUE(it.More());
  const TopoDS_Face   bridge = TopoDS::Face(it.Current());
  BRepAdaptor_Surface surface(bridge);
  const double        u0 = surface.FirstUParameter(), u1 = surface.LastUParameter(),
               v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
  // The blend runs along u or v; find the direction in which it reaches the wall by its ends.
  const bool alongU = surface.Value(u0, 0.5 * (v0 + v1)).Distance(gp_Pnt(10, 3, 0)) < 1e-6
                      || surface.Value(u1, 0.5 * (v0 + v1)).Distance(gp_Pnt(10, 3, 0)) < 1e-6;
  const auto at = [&](double t, double s) {
    return alongU ? surface.Value(u0 + (u1 - u0) * t, v0 + (v1 - v0) * s)
                  : surface.Value(u0 + (u1 - u0) * s, v0 + (v1 - v0) * t);
  };
  const auto normal = [&](double t, double s) {
    return alongU ? normalAt(bridge, u0 + (u1 - u0) * t, v0 + (v1 - v0) * s)
                  : normalAt(bridge, u0 + (u1 - u0) * s, v0 + (v1 - v0) * t);
  };
  for (const double s : {0.0, 0.25, 0.5, 1.0})
  {
    const bool startsOnFloor = at(0, s).Distance(gp_Pnt(10, 6 * s, 0)) < 1e-6
                               || at(0, s).Distance(gp_Pnt(10, 6 * (1 - s), 0)) < 1e-6;
    const gp_Pnt start  = startsOnFloor ? at(0, s) : at(1, s),
                 end    = startsOnFloor ? at(1, s) : at(0, s);
    const gp_Vec nStart = normal(startsOnFloor ? 0 : 1, s), nEnd = normal(startsOnFloor ? 1 : 0, s);
    EXPECT_NEAR(start.Z(), 0, 1e-6);
    EXPECT_NEAR(start.X(), 10, 1e-6);
    EXPECT_NEAR(end.X(), 16, 1e-6);
    EXPECT_NEAR(end.Z(), 3, 1e-6);
    EXPECT_NEAR(std::abs(nStart.Z()), 1, 1e-7);
    EXPECT_NEAR(std::abs(nEnd.X()), 1, 1e-7);
  }
  for (const double t : {0.1, 0.4, 0.7})
  {
    const gp_Pnt p = at(t, 0.5);
    EXPECT_GE(p.X(), 10 - 1e-9);
    EXPECT_LE(p.X(), 16 + 1e-9);
    EXPECT_GE(p.Z(), -1e-9);
    EXPECT_LE(p.Z(), 3 + 1e-9);
  }
}

// A chamfer-shaped bridge (shape 1) is the ruled strip between the two edges: here between
// (10, y, 0) and (16, y, 3) for y in 0..6, a plane of area 6 sqrt(36 + 9).
TEST(BRepSurfaceExt, ChamferBridgeIsTheRuledStripBetweenTheEdges)
{
  const TopoDS_Face floor =
    planeFace(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0), 0, 10, 0, 6);
  const TopoDS_Face wall =
    planeFace(gp_Pnt(16, 0, 3), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0), 0, 6, 0, 10);
  const auto floorEdge =
    edgesWhere(floor, [](const double* l) { return near(l[0], 10) && near(l[3], 10); });
  const auto wallEdge =
    edgesWhere(wall, [](const double* l) { return near(l[2], 3) && near(l[5], 3); });
  const auto          result = BRepExt::bridgeSurface(session(),
                                             BRepExt::store(session(), floorEdge[0]),
                                             BRepExt::store(session(), floor),
                                             BRepExt::store(session(), wallEdge[0]),
                                             BRepExt::store(session(), wall),
                                             1,
                                             1,
                                             1.0,
                                             1.0,
                                             1e-7,
                                             1);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  EXPECT_NEAR(area(made), 6 * std::sqrt(45.0), 1e-6 * 6 * std::sqrt(45.0));
}

// Aligning a tilted patch to the floor along their shared edge makes the patch tangent to the
// floor there. The patch is a bicubic B-spline over x in 10..16, y in 0..6 whose height grows like
// 0.5 i + 0.1 i j along x, so it meets the floor z = 0 at an angle at x = 10 (normal tilted by
// atan 0.5). After a G1 alignment the normals along the edge are (0, 0, 1) to 1e-9 and the edge
// itself has not moved: alignment must change the second row of poles only.
TEST(BRepSurfaceExt, G1AlignmentMakesAPatchTangentToTheFloorAndKeepsTheEdge)
{
  NCollection_Array2<gp_Pnt> poles(1, 4, 1, 4);
  for (int i = 1; i <= 4; ++i)
    for (int j = 1; j <= 4; ++j)
      poles(i, j) =
        gp_Pnt(10 + 2.0 * (i - 1), 2.0 * (j - 1), 0.5 * (i - 1) + 0.1 * (i - 1) * (j - 1));
  NCollection_Array1<double> knots(1, 2);
  NCollection_Array1<int>    mult(1, 2);
  knots(1) = 0;
  knots(2) = 1;
  mult(1)  = 4;
  mult(2)  = 4;
  const TopoDS_Face patch =
    BRepBuilderAPI_MakeFace(new Geom_BSplineSurface(poles, knots, knots, mult, mult, 3, 3), 1e-7)
      .Face();
  const TopoDS_Face floor =
    planeFace(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0), 0, 10, 0, 6);
  const auto patchEdge =
    edgesWhere(patch, [](const double* l) { return near(l[0], 10) && near(l[3], 10); });
  const auto floorEdge =
    edgesWhere(floor, [](const double* l) { return near(l[0], 10) && near(l[3], 10); });
  ASSERT_EQ(patchEdge.size(), 1u);
  ASSERT_EQ(floorEdge.size(), 1u);
  const auto          result = BRepExt::alignSurface(session(),
                                            BRepExt::store(session(), patch),
                                            BRepExt::store(session(), patchEdge[0]),
                                            BRepExt::store(session(), floorEdge[0]),
                                            BRepExt::store(session(), floor),
                                            1,
                                            1e-7,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            false,
                                            0,
                                            1,
                                            false,
                                            0,
                                            0);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  TopExp_Explorer     it(made, TopAbs_FACE);
  ASSERT_TRUE(it.More());
  const TopoDS_Face   aligned = TopoDS::Face(it.Current());
  BRepAdaptor_Surface surface(aligned);
  for (const double t : {0.0, 0.3, 0.6, 1.0})
  {
    // The edge at x = 10 is u = 0 (or u = 1 if the face was turned); find it by position.
    const double u0 = surface.FirstUParameter(), u1 = surface.LastUParameter(),
                 v0 = surface.FirstVParameter(), v1 = surface.LastVParameter();
    const bool uEdge = std::abs(surface.Value(u0, v0).X() - 10) < 1e-6
                       && std::abs(surface.Value(u0, v1).X() - 10) < 1e-6;
    const double u = uEdge ? u0 : u0 + (u1 - u0) * t, v = uEdge ? v0 + (v1 - v0) * t : v0;
    const gp_Pnt p = surface.Value(u, v);
    EXPECT_NEAR(p.X(), 10, 1e-9);
    EXPECT_NEAR(p.Z(), 0, 1e-9);
    EXPECT_NEAR(std::abs(normalAt(aligned, u, v).Z()), 1, 1e-9);
  }
}

// Two edges that meet, a line of length 10 along x and a quarter arc of radius 5 in the plane x = 0
// from the same point, are completed to a four-sided surface by translating each along the other.
// That surface is the cylinder patch swept by the arc along the line: its area is 10 times the
// arc's length 5 pi / 2, that is 25 pi, and every point is at
// distance 5 from the axis parallel to x through (0, 0, 5). A completion that closes the wrong
// corner, or a surface that bulges between the boundaries, changes the area or the radius.
TEST(BRepSurfaceExt, SquareSurfaceOfALineAndAnArcIsTheCylinderPatchTheyBound)
{
  const double  r = 5, tolerance = 1e-4;
  const gp_Ax2  frame(gp_Pnt(0, 0, r), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0));
  const gp_Circ circle(frame, r);
  // From the lowest point (0, 0, 0) a quarter turn to (0, r, r).
  const TopoDS_Edge arc =
    BRepBuilderAPI_MakeEdge(
      GC_MakeArcOfCircle(circle, gp_Pnt(0, 0, 0), gp_Pnt(0, r, r), true).Value())
      .Edge();
  const TopoDS_Edge line = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(10, 0, 0)).Edge();
  const auto        result =
    BRepExt::squareSurface(session(),
                           {BRepExt::store(session(), line), BRepExt::store(session(), arc)},
                           {0, 0},
                           {0, 0},
                           3,
                           5,
                           1,
                           2,
                           tolerance);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_NEAR(area(made), 25 * M_PI, 1e-3 * 25 * M_PI);
  double worst = 0;
  for (TopExp_Explorer it(made, TopAbs_FACE); it.More(); it.Next())
  {
    BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; j <= 8; ++j)
      {
        const gp_Pnt p =
          surface.Value(surface.FirstUParameter()
                          + (surface.LastUParameter() - surface.FirstUParameter()) * i / 8,
                        surface.FirstVParameter()
                          + (surface.LastVParameter() - surface.FirstVParameter()) * j / 8);
        worst = std::max(worst, std::abs(std::hypot(p.Y(), p.Z() - r) - r));
      }
  }
  EXPECT_LT(worst, tolerance) << "radial error " << worst;
}

namespace
{

//! A bicubic B-spline face of 4 x 4 spans over a square of side `side`, its height a smooth wave
//! of amplitude `wave`, placed by `placement`.
TopoDS_Face waveFace(double side, double wave, const gp_Trsf& placement)
{
  const int                  spans = 4, count = spans + 3;
  NCollection_Array2<gp_Pnt> poles(1, count, 1, count);
  for (int i = 1; i <= count; ++i)
    for (int j = 1; j <= count; ++j)
      poles(i, j) = gp_Pnt(side * (i - 1) / (count - 1),
                           side * (j - 1) / (count - 1),
                           wave * std::sin(1.3 * i) * std::cos(0.9 * j));
  NCollection_Array1<double> knots(1, spans + 1);
  NCollection_Array1<int>    mult(1, spans + 1);
  for (int k = 1; k <= spans + 1; ++k)
  {
    knots(k) = double(k - 1) / spans;
    mult(k)  = (k == 1 || k == spans + 1) ? 4 : 1;
  }
  Handle(Geom_BSplineSurface) surface =
    new Geom_BSplineSurface(poles, knots, knots, mult, mult, 3, 3);
  surface->Transform(placement);
  return BRepBuilderAPI_MakeFace(surface, 1e-7).Face();
}

} // namespace

// The face separation pre-pass must never prove two faces clear of each other that are not: a
// pair it reports as proven has a positive minimum distance, measured by BRepExtrema. Two wavy
// 10 x 10 patches are put in 40 poses (a grid of offsets in x, z and a tilt) from far apart
// through touching to crossing; every pair declared proven is checked against the distance, and
// the poses include both proven and unproven pairs, so that neither "proves everything" nor
// "proves nothing" passes.
TEST(BRepSurfaceExt, SeparationNeverProvesFacesThatTouchOrCross)
{
  const TopoDS_Face base   = waveFace(10, 1.5, gp_Trsf());
  int               proven = 0, unproven = 0;
  for (int pose = 0; pose < 40; ++pose)
  {
    gp_Trsf tilt, shift;
    tilt.SetRotation(gp_Ax1(gp_Pnt(5, 5, 0), gp_Dir(1, 0.3, 0)), 0.15 * (pose % 5));
    shift.SetTranslation(gp_Vec(2.5 * (pose % 8) - 5, 0.7 * (pose % 3), 0.9 * (pose / 8) - 1.8));
    const TopoDS_Face other = waveFace(10, 1.5, shift * tilt);
    BRep_Builder      builder;
    TopoDS_Compound   both;
    builder.MakeCompound(both);
    builder.Add(both, base);
    builder.Add(both, other);
    BRepExt::SeparationCount count;
    BRepExt::unprovenParts(both, &count);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> faces;
    TopExp::MapShapes(both, TopAbs_FACE, faces);
    ASSERT_EQ(faces.Extent(), 2);
    bool crossProven = false;
    for (const auto& pair : count.provenPairs)
      crossProven = crossProven || pair[0] != pair[1];
    if (crossProven)
    {
      ++proven;
      BRepExtrema_DistShapeShape distance(base, other);
      ASSERT_TRUE(distance.IsDone());
      EXPECT_GT(distance.Value(), 1e-9)
        << "pose " << pose << " proven but at distance " << distance.Value();
    }
    else
      ++unproven;
  }
  EXPECT_GT(proven, 0);
  EXPECT_GT(unproven, 0);
}

// A non-planar quintic/cubic net beside a plane. Only the first continuity rows
// and one extra row may move; compare its regularisation path at several poses/scales.
TEST(BRepSurfaceExt, AlignmentPolePriorPreservesTheInputWithoutRelaxingTheSeam)
{
  for (const double scale : {0.01, 1.0, 100.0})
    for (const bool posed : {false, true})
      for (const int order : {1, 2})
      {
        gp_Trsf pose;
        if (posed)
        {
          pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.7);
          pose.SetTranslationPart(gp_Vec(7 * scale, -3 * scale, 5 * scale));
        }
        const auto point = [&](double x, double y, double z) {
          return gp_Pnt(scale * x, scale * y, scale * z).Transformed(pose);
        };
        NCollection_Array2<gp_Pnt> poles(1, 6, 1, 4);
        for (int i = 1; i <= 6; ++i)
          for (int j = 1; j <= 4; ++j)
            poles(i, j) = point(10 + 2 * (i - 1), 2 * (j - 1),
                               0.3 * (i - 1) + 0.04 * (i - 1) * (i - 1)
                                 + 0.02 * (i - 1) * (j - 1));
        NCollection_Array1<double> knots(1, 2);
        NCollection_Array1<int> mu(1, 2), mv(1, 2);
        knots(1) = 0; knots(2) = 1;
        mu(1) = mu(2) = 6;
        mv(1) = mv(2) = 4;
        const occ::handle<Geom_BSplineSurface> input =
          new Geom_BSplineSurface(poles, knots, knots, mu, mv, 5, 3);
        const TopoDS_Face patch = BRepBuilderAPI_MakeFace(input, 1e-8 * scale).Face();
        const gp_Dir normal = gp_Dir(0, 0, 1).Transformed(pose);
        const TopoDS_Face floor = planeFace(point(0, 0, 0), normal,
                                           gp_Dir(1, 0, 0).Transformed(pose),
                                           0, 10 * scale, 0, 6 * scale);
        const auto seam = [&](const TopoDS_Face& face) {
          for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
          {
            BRepAdaptor_Curve c(TopoDS::Edge(it.Current()));
            const gp_Pnt mid = c.Value((c.FirstParameter() + c.LastParameter()) / 2);
            if (mid.Distance(point(10, 3, 0)) < 1e-6 * scale)
              return TopoDS::Edge(it.Current());
          }
          return TopoDS_Edge();
        };
        const auto faceId = BRepExt::store(session(), patch),
                   edgeId = BRepExt::store(session(), seam(patch)),
                   refId = BRepExt::store(session(), seam(floor)),
                   supportId = BRepExt::store(session(), floor);
        const auto run = [&](int extra, double w) {
          const auto result = BRepExt::alignSurface(session(), faceId, edgeId, refId, supportId,
            order, 1e-7 * scale, extra, 0, 0, 0, 0, false, 0, 1, false, 0, 0, w);
          const auto face = TopoDS::Face(BRepExt::get(session(), result.resultId));
          EXPECT_TRUE(BRepExt::invalidity(face).empty());
          return occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(face));
        };
        const auto baseline = run(0, 0), inert = run(0, 1);
        for (int i = 1; i <= 6; ++i)
          for (int j = 1; j <= 4; ++j)
            EXPECT_TRUE(baseline->Pole(i, j).IsEqual(inert->Pole(i, j), 0));
        double previous = std::numeric_limits<double>::infinity();
        double first = 0, last = 0;
        for (const double w : {0.0, 0.1, 0.5, 1.0})
        {
          const auto made = run(1, w);
          double movement = 0;
          for (int i = 1; i <= 6; ++i)
            for (int j = 1; j <= 4; ++j)
            {
              movement += made->Pole(i, j).SquareDistance(poles(i, j));
              EXPECT_TRUE(input->Pole(i, j).IsEqual(poles(i, j), 0));
              if (i > order + 2)
                EXPECT_TRUE(made->Pole(i, j).IsEqual(poles(i, j), 0));
            }
          EXPECT_LE(movement, previous + 1e-9 * scale * scale);
          previous = movement;
          if (w == 0) first = movement;
          last = movement;
          for (int j = 0; j <= 32; ++j)
          {
            gp_Pnt p;
            gp_Vec du, dv, duu, dvv, duv;
            made->D2(0, j / 32.0, p, du, dv, duu, dvv, duv);
            EXPECT_LT(std::abs(gp_Vec(point(0, 0, 0), p).Dot(gp_Vec(normal))), 1e-6 * scale);
            EXPECT_LT(du.Crossed(dv).Normalized().Crossed(gp_Vec(normal)).Magnitude(), 1e-6);
            if (order == 2)
              EXPECT_LT(std::abs(duu.Dot(gp_Vec(normal))), 1e-5 * scale);
          }
        }
        EXPECT_LT(last, first * 0.999);
        EXPECT_THROW(run(1, -0.1), BRepExt::Failure);
        EXPECT_THROW(run(1, std::numeric_limits<double>::quiet_NaN()), BRepExt::Failure);
        for (int i = 1; i <= 6; ++i)
          for (int j = 1; j <= 4; ++j)
            EXPECT_TRUE(input->Pole(i, j).IsEqual(poles(i, j), 0));
      }
}

namespace
{
struct CollinearPair
{
  TopoDS_Face moving, reference;
  TopoDS_Edge movingEdge, referenceEdge;
};

CollinearPair collinearPair(double scale,
                            bool   posed,
                            int    degree,
                            bool   rational,
                            bool   folded = false)
{
  gp_Trsf pose;
  if (posed)
  {
    pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.7);
    pose.SetTranslationPart(gp_Vec(2 * scale, -3 * scale, 4 * scale));
  }
  const auto make = [&](bool moving) {
    const int                  count = degree == 1 ? 2 : 4;
    NCollection_Array2<gp_Pnt> poles(1, count, 1, count);
    NCollection_Array2<double> weights(1, count, 1, count);
    NCollection_Array1<double> knots(1, 2);
    NCollection_Array1<int>    mult(1, 2);
    knots(1) = 0;
    knots(2) = 1;
    mult(1) = mult(2) = count;
    for (int i = 1; i <= count; ++i)
      for (int j = 1; j <= count; ++j)
      {
        const double x = (moving && !folded ? 1 : -1) * (i - 1.0);
        poles(i, j) =
          gp_Pnt(
            scale * x,
            scale * (j - 1),
            scale * (degree == 1 ? 0 : 0.1 * (j - 1) * (j - 1) + (moving && !folded ? 0.2 * x : 0)))
            .Transformed(pose);
        weights(i, j) = rational ? 1 + (moving ? 0.11 : 0.07) * j : 1;
      }
    occ::handle<Geom_BSplineSurface> s =
      new Geom_BSplineSurface(poles, weights, knots, knots, mult, mult, count - 1, count - 1);
    s->IncreaseDegree(degree, !moving && degree == 6 ? 9 : degree);
    if (degree > 1)
    {
      if (moving)
        s->InsertVKnot(0.37, 1, 1e-12);
      else
        s->InsertVKnot(0.61, 2, 1e-12);
    }
    if (posed)
      s->UReverse();
    if (!moving && posed)
      s->VReverse();
    if (degree == 9)
      s->ExchangeUV();
    auto face = BRepBuilderAPI_MakeFace(s, 1e-7).Face();
    // The same relative positional allowance at every homothetic scale.
    const double allowance = std::max(1e-7, 1e-7 * scale);
    BRep_Builder topology;
    topology.UpdateFace(face, allowance);
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
      topology.UpdateEdge(TopoDS::Edge(it.Current()), allowance);
    for (TopExp_Explorer it(face, TopAbs_VERTEX); it.More(); it.Next())
      topology.UpdateVertex(TopoDS::Vertex(it.Current()), allowance);
    return face;
  };
  CollinearPair pair;
  pair.moving     = make(true);
  pair.reference  = make(false);
  const auto edge = [&](const TopoDS_Face& f) {
    for (TopExp_Explorer it(f, TopAbs_EDGE); it.More(); it.Next())
      if (BRepExt::surfacing::boundarySide("test", f, TopoDS::Edge(it.Current()))
          == (degree == 9 ? 2 : 0) + (posed ? 1 : 0))
        return TopoDS::Edge(it.Current());
    return TopoDS_Edge();
  };
  pair.movingEdge    = edge(pair.moving);
  pair.referenceEdge = edge(pair.reference);
  return pair;
}

void checkCollinear(double scale, bool posed, int degree, bool rational, bool symmetric)
{
  SCOPED_TRACE(::testing::Message()
               << "scale=" << scale << " posed=" << posed << " degree=" << degree
               << " rational=" << rational << " symmetric=" << symmetric);
  const auto pair = collinearPair(scale, posed, degree, rational);
  ASSERT_TRUE(BRepExt::invalidity(pair.moving).empty());
  ASSERT_TRUE(BRepExt::invalidity(pair.reference).empty());
  const auto result = BRepExt::alignCollinear(session(),
                                              BRepExt::store(session(), pair.moving),
                                              BRepExt::store(session(), pair.movingEdge),
                                              BRepExt::store(session(), pair.referenceEdge),
                                              BRepExt::store(session(), pair.reference),
                                              symmetric,
                                              true,
                                              0,
                                              1,
                                              0,
                                              0,
                                              1e-7 * scale);
  std::vector<TopoDS_Face> faces;
  for (TopExp_Explorer it(BRepExt::get(session(), result.resultId), TopAbs_FACE); it.More();
       it.Next())
    faces.push_back(TopoDS::Face(it.Current()));
  ASSERT_EQ(faces.size(), symmetric ? 2u : 1u);
  for (const auto& face : faces)
  {
    EXPECT_TRUE(BRepExt::invalidity(face).empty());
    const double allowance = std::max(1e-7, 1e-7 * scale);
    EXPECT_LE(BRep_Tool::Tolerance(face), allowance);
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
      EXPECT_LE(BRep_Tool::Tolerance(TopoDS::Edge(it.Current())), allowance);
    for (TopExp_Explorer it(face, TopAbs_VERTEX); it.More(); it.Next())
      EXPECT_LE(BRep_Tool::Tolerance(TopoDS::Vertex(it.Current())), allowance);
  }
  const auto a = BRep_Tool::Surface(faces[0]);
  const auto b = BRep_Tool::Surface(symmetric ? faces[1] : pair.reference);
  for (int j = 0; j <= 40; ++j)
  {
    gp_Pnt       pa, pb;
    gp_Vec       au, av, bu, bv;
    const double across = posed ? 1 : 0;
    const double t = j / 40.0, rt = posed ? 1 - t : t;
    a->D1(degree == 9 ? t : across, degree == 9 ? across : t, pa, au, av);
    b->D1(degree == 9 ? rt : across, degree == 9 ? across : rt, pb, bu, bv);
    const auto   input = BRep_Tool::Surface(pair.moving);
    const double far   = 1 - across;
    EXPECT_LE(a->Value(degree == 9 ? t : far, degree == 9 ? far : t)
                .Distance(input->Value(degree == 9 ? t : far, degree == 9 ? far : t)),
              scale * 1e-12);
    EXPECT_LE(pa.Distance(pb), scale * 1e-12);
    EXPECT_LE(au.Crossed(av).Normalized().Crossed(bu.Crossed(bv).Normalized()).Magnitude(), 1e-11);
  }
}
} // namespace

TEST(BRepSurfaceExt, CollinearSymmetricUsesAnAngularFoldBoundAtEveryScale)
{
  for (double scale : {0.001, 1.0, 3e9})
    for (bool posed : {false, true})
      for (int degree : {1, 3})
        EXPECT_NO_THROW(checkCollinear(scale, posed, degree, false, true))
          << "scale=" << scale << " posed=" << posed << " degree=" << degree;
}

TEST(BRepSurfaceExt, CollinearHarmonisesHigherDegreesAndRationalRowsExactly)
{
  for (double scale : {0.01, 1.0, 100.0})
    for (bool posed : {false, true})
      for (int degree : {6, 9})
        for (bool rational : {false, true})
          for (bool symmetric : {false, true})
            EXPECT_NO_THROW(checkCollinear(scale, posed, degree, rational, symmetric));
}

TEST(BRepSurfaceExt, CollinearConstantRatioChecksRationalDerivativeWeights)
{
  const auto pair = collinearPair(1, false, 6, true);
  auto       s    = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(pair.reference)->Copy());
  s->SetWeight(2, 2, 2 * s->Weight(2, 2));
  const auto  reference = BRepBuilderAPI_MakeFace(s, 1e-7).Face();
  TopoDS_Edge edge;
  for (TopExp_Explorer it(reference, TopAbs_EDGE); it.More(); it.Next())
    if (BRepExt::surfacing::boundarySide("test", reference, TopoDS::Edge(it.Current())) == 0)
      edge = TopoDS::Edge(it.Current());
  try
  {
    BRepExt::alignCollinear(session(),
                            BRepExt::store(session(), pair.moving),
                            BRepExt::store(session(), pair.movingEdge),
                            BRepExt::store(session(), edge),
                            BRepExt::store(session(), reference),
                            false,
                            true,
                            0,
                            1,
                            0,
                            0,
                            1e-7);
    FAIL() << "a Euclidean row ratio does not certify general rational derivatives";
  }
  catch (const BRepExt::Failure& failure)
  {
    EXPECT_EQ(failure.reason, BRepExt::unsupported_geometry);
  }
  // General positive rational weights still admit the pole-collinearity option.
  // Only its boundary equality is promised when the ratio is allowed to vary.
  const auto result = BRepExt::alignCollinear(session(),
                                              BRepExt::store(session(), pair.moving),
                                              BRepExt::store(session(), pair.movingEdge),
                                              BRepExt::store(session(), edge),
                                              BRepExt::store(session(), reference),
                                              false,
                                              false,
                                              0,
                                              1,
                                              0,
                                              0,
                                              1e-7);
  const auto made   = TopoDS::Face(BRepExt::get(session(), result.resultId));
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const auto surface = BRep_Tool::Surface(made);
  for (int j = 0; j <= 40; ++j)
    EXPECT_LE(surface->Value(0, j / 40.0).Distance(s->Value(0, j / 40.0)), 1e-12);
}

TEST(BRepSurfaceExt, CollinearStillRejectsOpposingRowsAtEveryScale)
{
  for (double scale : {0.001, 1.0, 3e9})
    for (bool posed : {false, true})
    {
      const auto pair = collinearPair(scale, posed, 3, false, true);
      EXPECT_THROW(BRepExt::alignCollinear(session(),
                                           BRepExt::store(session(), pair.moving),
                                           BRepExt::store(session(), pair.movingEdge),
                                           BRepExt::store(session(), pair.referenceEdge),
                                           BRepExt::store(session(), pair.reference),
                                           true,
                                           true,
                                           0,
                                           1,
                                           0,
                                           0,
                                           1e-7),
                   BRepExt::Failure);
    }
}
