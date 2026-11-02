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
// shared boundary.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
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
