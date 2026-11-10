// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepAlgoExt_TrimDomain.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepAlgoExt_ReadBack.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gtest/gtest.h>
#include <cmath>

// Exact moments of a rectangle with an eccentric circular hole. The independent
// signed subtraction catches ignoring inner wires, wrong winding and UV Jacobian.
TEST(BRepAlgoExtShared, EccentricHoleMoments)
{
  for (double s : {.1, 1., 10.})
  {
    BRepBuilderAPI_MakeFace maker(gp_Pln(gp_Pnt(0, 0, 2 * s), gp::DZ()), 0, 4 * s, 0, 3 * s);
    const auto              hole =
      BRepBuilderAPI_MakeWire(
        BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(s, s, 2 * s), gp::DZ()), .4 * s)).Edge())
        .Wire();
    maker.Add(TopoDS::Wire(hole.Reversed()));
    const gp_Pnt origin(0, 0, 0);
    const auto   result = BRepExt::trimmedProperties(maker.Face(), &origin);
    const double disk = M_PI * .16 * s * s, area = 12 * s * s - disk;
    ASSERT_TRUE(result.converged);
    EXPECT_NEAR(result.area, area, 1e-9 * s * s);
    EXPECT_NEAR(result.center.X(), (24 * s * s * s - disk * s) / area, 1e-9 * s);
    EXPECT_NEAR(result.center.Y(), (18 * s * s * s - disk * s) / area, 1e-9 * s);
    EXPECT_NEAR(result.center.Z(), 2 * s, 1e-9 * s);
    EXPECT_NEAR(result.volume, area * 2 * s / 3, 1e-9 * s * s * s);
  }
}

// Stored tolerance is a maximum over topology, without a document floor.
// A nested vertex, not the root container or its first child, holds the maximum.
TEST(BRepAlgoExtShared, NestedToleranceMaximum)
{
  BRep_Builder    b;
  TopoDS_Compound outer, inner;
  b.MakeCompound(outer);
  b.MakeCompound(inner);
  for (double t : {1e-5, 3e-4, 2e-5})
  {
    TopoDS_Vertex v;
    b.MakeVertex(v, gp_Pnt(t, 0, 0), t);
    b.Add(inner, v);
  }
  b.Add(outer, inner);
  EXPECT_DOUBLE_EQ(BRepExt::tolerance::own(outer), 0);
  EXPECT_DOUBLE_EQ(BRepExt::tolerance::largest(outer), 3e-4);
}

// The reader normalizes coordinates again. These non-axis-aligned direction
// families catch returning an unstable input or perturbing the direction grossly.
TEST(BRepAlgoExtShared, DirectionReaderFixedPoints)
{
  for (int i = 1; i <= 128; ++i)
  {
    const gp_Dir input(std::sin(i * .37), std::cos(i * .19), std::sin(i * .11) + .4);
    const auto   settled = BRepExt::readBack(input);
    gp_Dir       read;
    read.SetCoord(settled.X(), settled.Y(), settled.Z());
    EXPECT_TRUE(BRepExt::same(settled, read));
    EXPECT_LT(gp_Vec(input).Crossed(gp_Vec(settled)).Magnitude(), 1e-12);
  }
}

#include <BRepAlgoExt_Properties.hxx>
#include <gp_Sphere.hxx>

// Integral r^2 cos(v) du dv on a sphere rectangle; probes both closed-form
// and adaptive paths without deriving expected values from another OCCT measure.
TEST(BRepAlgoExtShared, SphereRectangleArea)
{
  for (double r : {.1, 1., 10.})
  {
    const auto   face = BRepBuilderAPI_MakeFace(gp_Sphere(gp_Ax3(), r), .2, 1.3, -.4, .7).Face();
    const double expected = r * r * 1.1 * (std::sin(.7) - std::sin(-.4));
    EXPECT_NEAR(BRepExt::faceArea(face), expected, 1e-10 * r * r);
    bool converged = true;
    EXPECT_NEAR(BRepExt::faceMeasure(face, true, gp_Pnt(), &converged), expected, 1e-8 * r * r);
    EXPECT_TRUE(converged);
  }
}

#include <Geom_BezierCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>

// C=(u,u^2,1), D=(0,0,1) gives flux u^2/3. The curved upper
// trim v=1+u^2 gives integral (1/3+1/5)/3=8/45. Shortening that
// trim leaves a tolerance-sized endpoint gap; its closing chord changes
// the integral only cubically in the gap, below binary64 rounding here.
TEST(BRepAlgoExtShared, PolynomialExtrusionFluxWithCurvedTrim)
{
  for (bool gap : {false, true})
  {
    NCollection_Array1<gp_Pnt> poles(1, 3);
    poles(1) = gp_Pnt(0, 0, 1);
    poles(2) = gp_Pnt(.5, 0, 1);
    poles(3) = gp_Pnt(1, 1, 1);
    const Handle(Geom_SurfaceOfLinearExtrusion) surface =
      new Geom_SurfaceOfLinearExtrusion(new Geom_BezierCurve(poles), gp::DZ());
    auto line = [&](const gp_Pnt2d& a, const gp_Pnt2d& b) {
      return BRepBuilderAPI_MakeEdge(new Geom2d_Line(a, gp_Dir2d(gp_Vec2d(a, b))),
                                     surface,
                                     0,
                                     a.Distance(b))
        .Edge();
    };
    NCollection_Array1<gp_Pnt2d> top(1, 3);
    top(1)            = gp_Pnt2d(1, 2);
    top(2)            = gp_Pnt2d(.5, 1);
    top(3)            = gp_Pnt2d(0, 1);
    const auto curved = BRepBuilderAPI_MakeEdge(new Geom2d_BezierCurve(top), surface, 0, 1).Edge();
    BRepBuilderAPI_MakeWire wire;
    wire.Add(line(gp_Pnt2d(0, 0), gp_Pnt2d(1, 0)));
    wire.Add(line(gp_Pnt2d(1, 0), gp_Pnt2d(1, 2)));
    wire.Add(curved);
    wire.Add(line(gp_Pnt2d(0, 1), gp_Pnt2d(0, 0)));
    const auto face = BRepBuilderAPI_MakeFace(surface, wire.Wire(), true).Face();
    if (gap)
      BRep_Builder().Range(curved, face, 8.6e-10, 1);
    for (double z : {-100., 0., 100.})
    {
      bool converged = true;
      EXPECT_NEAR(BRepExt::faceMeasure(face, false, gp_Pnt(0, 0, z), &converged), 8. / 45, 3e-16);
      EXPECT_TRUE(converged);
      EXPECT_NEAR(
        BRepExt::faceMeasure(TopoDS::Face(face.Reversed()), false, gp_Pnt(0, 0, z), nullptr),
        -8. / 45,
        3e-16);
    }
  }
}

// With degree-nine Bezier poles X_8=a, X_9=Y_9=1 and all others
// zero, det(C,C',D)=9a*u^16. Its rectangle flux is exactly 3a/17.
// Eight Gauss points cannot integrate degree 16 exactly; nine can.
TEST(BRepAlgoExtShared, HighDegreePolynomialExtrusionRectangleFlux)
{
  NCollection_Array1<gp_Pnt> poles(1, 10);
  for (int i = 1; i <= 10; ++i)
    poles(i) = gp_Pnt(0, 0, 1);
  const double a = 1. / 9;
  poles(9)       = gp_Pnt(a, 0, 1);
  poles(10)      = gp_Pnt(1, 1, 1);
  const Handle(Geom_SurfaceOfLinearExtrusion) surface =
    new Geom_SurfaceOfLinearExtrusion(new Geom_BezierCurve(poles), gp::DZ());
  const auto face      = BRepBuilderAPI_MakeFace(surface, 0, 1, 0, 1, 1e-7).Face();
  bool       converged = true;
  EXPECT_NEAR(BRepExt::faceMeasure(face, false, gp_Pnt(), &converged), 3 * a / 17, 2e-17);
  EXPECT_TRUE(converged);
}

#include <Geom_BSplineCurve.hxx>
#include <BRepLib.hxx>
#include <BRepCheck_Analyzer.hxx>

// The periodic square carrier has g=1/3 on two of its four spans and
// zero on the others. The roof v=1+(u-shift)/4 integrates to one on
// these seam-crossing charts. A shortened linear trim has the same closure.
TEST(BRepAlgoExtShared, PeriodicPolynomialExtrusionFlux)
{
  NCollection_Array1<gp_Pnt> poles(1, 4);
  poles(1) = gp_Pnt(0, 0, 1);
  poles(2) = gp_Pnt(1, 0, 1);
  poles(3) = gp_Pnt(1, 1, 1);
  poles(4) = gp_Pnt(0, 1, 1);
  NCollection_Array1<double> knots(1, 5);
  NCollection_Array1<int>    mults(1, 5);
  for (int i = 1; i <= 5; ++i)
  {
    knots(i) = i - 1;
    mults(i) = 1;
  }
  const Handle(Geom_SurfaceOfLinearExtrusion) surface =
    new Geom_SurfaceOfLinearExtrusion(new Geom_BSplineCurve(poles, knots, mults, 1, true),
                                      gp::DZ());
  for (int shift : {0, 2, -2})
    for (bool gap : {false, true})
    {
      const std::array<gp_Pnt2d, 4> uv{gp_Pnt2d(shift, 0),
                                       gp_Pnt2d(shift + 4, 0),
                                       gp_Pnt2d(shift + 4, 2),
                                       gp_Pnt2d(shift, 1)};
      BRepBuilderAPI_MakeWire       wire;
      TopoDS_Edge                   top;
      for (int i = 0; i < 4; ++i)
      {
        NCollection_Array1<gp_Pnt2d> p(1, 2);
        p(1)            = uv[i];
        p(2)            = uv[(i + 1) % 4];
        const auto edge = BRepBuilderAPI_MakeEdge(new Geom2d_BezierCurve(p), surface, 0, 1).Edge();
        wire.Add(edge);
        if (i == 2)
          top = edge;
      }
      const auto face = BRepBuilderAPI_MakeFace(surface, wire.Wire(), true).Face();
      BRepLib::BuildCurves3d(face);
      if (gap)
        BRep_Builder().Range(top, face, 8.6e-10, 1);
      ASSERT_TRUE(BRepCheck_Analyzer(face).IsValid());
      bool converged = true;
      EXPECT_DOUBLE_EQ(BRepExt::faceMeasure(face, false, gp_Pnt(), &converged), 1);
      EXPECT_TRUE(converged);
    }
}
