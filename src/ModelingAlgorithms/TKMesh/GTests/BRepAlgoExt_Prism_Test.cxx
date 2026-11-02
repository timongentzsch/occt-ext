// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepAlgoExt_Prism.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GeomConvert.hxx>
#include <Geom_Circle.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pln.hxx>
#include <gtest/gtest.h>

// The inradius of a circular disk is exactly its radius. This catches a missing
// factor of two or incorrect transformed cap bounds, over scales and placements.
TEST(BRepAlgoExtPrism, CircularCapInradius)
{
  for (double radius : {.1, 1., 10.})
  {
    const gp_Ax2              frame(gp_Pnt(3 * radius, -2 * radius, 7 * radius), gp_Dir(1, 2, 3));
    const Handle(Geom_Circle) circle = new Geom_Circle(frame, radius);
    const auto                spline = GeomConvert::CurveToBSplineCurve(circle);
    const auto                edge   = BRepBuilderAPI_MakeEdge(spline).Edge();
    const auto cap  = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(edge).Wire(), true).Face();
    const auto body = BRepPrimAPI_MakePrism(cap, gp_Vec(frame.Direction()) * (4 * radius)).Shape();
    EXPECT_FALSE(BRepExt::blendConsumesPrismCap(body, {edge}, .5 * radius));
    EXPECT_TRUE(BRepExt::blendConsumesPrismCap(body, {edge}, 1.5 * radius));
  }
}

// The inward vector of every rectangle edge points toward the centre, even if
// the face orientation is reversed. This catches orientation-dependent normals.
TEST(BRepAlgoExtPrism, InwardIsIndependentOfFaceOrientation)
{
  for (double scale : {.1, 1., 10.})
  {
    const auto face = BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()), 0, 2 * scale, 0, scale).Face();
    for (bool reverse : {false, true})
    {
      const auto support = reverse ? TopoDS::Face(face.Reversed()) : face;
      for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
      {
        const auto        edge = TopoDS::Edge(it.Current());
        BRepAdaptor_Curve curve(edge);
        const double      t = (curve.FirstParameter() + curve.LastParameter()) / 2;
        gp_Pnt            p;
        gp_Vec            tangent, inward;
        curve.D1(t, p, tangent);
        ASSERT_TRUE(BRepExt::supportInward(support, edge, t, tangent, inward));
        EXPECT_GT(inward.Dot(gp_Vec(p, gp_Pnt(scale, .5 * scale, 0))), 0);
      }
    }
  }
}
