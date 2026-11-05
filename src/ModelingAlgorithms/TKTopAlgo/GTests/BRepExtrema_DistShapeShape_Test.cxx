// Copyright (c) 2025 OPEN CASCADE SAS
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

#include <gtest/gtest.h>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepExtrema_ExtPC.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Trsf.hxx>

#include <cmath>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

//=================================================================================================

// Bug BUC60870: BRepExtrema_DistShapeShape should correctly compute the minimum distance
// between a linear edge and a vertex.
// The edge runs from (0,0,0) to (0,1,0); the vertex is at (0,0.3,1).
// The theoretical minimum distance is 1.0 (perpendicular from (0,0.3,1) to (0,0.3,0)).
TEST(BRepExtrema_DistShapeShapeTest, BUC60870_EdgeToVertexMinimumDistance)
{
  const TopoDS_Edge   anEdge  = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(0, 1, 0)).Edge();
  const TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(gp_Pnt(0, 0.3, 1)).Vertex();

  BRepExtrema_DistShapeShape aDist(anEdge, aVertex, 2.0);

  ASSERT_TRUE(aDist.IsDone()) << "BRepExtrema_DistShapeShape failed to compute";
  ASSERT_GE(aDist.NbSolution(), 1) << "Expected at least one solution";

  const double aTheorDist  = 1.0;
  const double aPercentMax = 0.01; // 1% relative tolerance
  EXPECT_NEAR(aDist.Value(), aTheorDist, aTheorDist * aPercentMax)
    << "Minimum distance deviates from expected value 1.0";
}

//=================================================================================================

// Verify that BRepExtrema_DistShapeShape does not crash when computing the distance
// between two edges and one of them has no 3D curve (null Geom_Curve from BRep_Tool::Curve).
// This exercises the null-check in the PERFORM_C0 helper.
TEST(BRepExtrema_DistShapeShapeTest, EdgeEdge_Null3DCurve_NoCrash)
{
  const TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(1, 0, 0)).Edge();

  const occ::handle<Geom_Circle> aCircle =
    new Geom_Circle(gp_Ax2(gp_Pnt(5, 0, 0), gp_Dir(0, 0, 1)), 1.0);
  TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(aCircle).Edge();

  BRep_Builder aBuilder;
  aBuilder.UpdateEdge(aEdge2, occ::handle<Geom_Curve>(), 1e-7);

  BRepExtrema_DistShapeShape aDist(aEdge1, aEdge2, 10.0);
  (void)aDist;
}

//=================================================================================================

// Verify that BRepExtrema_DistShapeShape does not crash when computing the distance
// between an edge with no 3D curve and a face.
// This exercises the null-check in Perform(Edge, Face).
TEST(BRepExtrema_DistShapeShapeTest, EdgeFace_Null3DCurve_NoCrash)
{
  const occ::handle<Geom_Circle> aCircle =
    new Geom_Circle(gp_Ax2(gp_Pnt(5, 0, 0), gp_Dir(0, 0, 1)), 1.0);
  TopoDS_Edge aEdge = BRepBuilderAPI_MakeEdge(aCircle).Edge();

  BRep_Builder aBuilder;
  aBuilder.UpdateEdge(aEdge, occ::handle<Geom_Curve>(), 1e-7);

  const TopoDS_Face aFace =
    BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 5), gp_Dir(0, 0, 1)), -10, 10, -10, 10);

  BRepExtrema_DistShapeShape aDist(aEdge, aFace, 10.0);
  (void)aDist;
}

// The normalized stationary-distance function has length units. A parameter
// tolerance must not make it "all null" on a short or slowly parameterized span.
TEST(BRepExtrema_DistShapeShapeTest, LinearSplineProjectionIndependentOfParameterScale)
{
  for (const double aLength : {1.0, 0.01, 0.001, 0.0001})
    for (const double aSpan : {1.0, 100.0, 10000.0})
      for (const double anAngle : {0.0, 0.73, 2.17})
        for (const bool isTrimmed : {false, true})
          for (const bool isRational : {false, true})
            {
              SCOPED_TRACE(::testing::Message() << aLength << " " << aSpan << " " << anAngle
                                               << " " << isTrimmed << " " << isRational);
              gp_Trsf aPose;
              aPose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), anAngle);
              aPose.SetTranslationPart(gp_Vec(2, -3, 5));
              NCollection_Array1<gp_Pnt> aPoles(1, 3);
              aPoles(1) = gp_Pnt(0, 0, 0);
              aPoles(2) = gp_Pnt(aLength, 0, 0);
              aPoles(3) = gp_Pnt(2 * aLength, aLength, 0);
              for (auto& aPole : aPoles)
                aPole.Transform(aPose);
              NCollection_Array1<double> aKnots(1, 3);
              aKnots(1) = 0;
              aKnots(2) = aSpan;
              aKnots(3) = 2 * aSpan;
              NCollection_Array1<int> aMults(1, 3);
              aMults(1) = aMults(3) = 2;
              aMults(2) = 1;
              NCollection_Array1<double> aWeights(1, 3);
              aWeights(1) = aWeights(3) = 1;
              aWeights(2) = 2;
              const occ::handle<Geom_BSplineCurve> aCurve = isRational
                ? new Geom_BSplineCurve(aPoles, aWeights, aKnots, aMults, 1)
                : new Geom_BSplineCurve(aPoles, aKnots, aMults, 1);
              const TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(
                aCurve, isTrimmed ? 0.1 * aSpan : 0.0, isTrimmed ? 0.8 * aSpan : 2 * aSpan);
              const gp_Pnt aQuery = gp_Pnt(0.37 * aLength, 1e-7, 0).Transformed(aPose);
              const TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(aQuery);
              BRepExtrema_DistShapeShape aDistance(aVertex, anEdge);
              ASSERT_TRUE(aDistance.IsDone());
              EXPECT_NEAR(aDistance.Value(), 1e-7, 4e-15);
              BRepExtrema_ExtPC anExtrema(aVertex, anEdge);
              ASSERT_TRUE(anExtrema.IsDone());
              ASSERT_GE(anExtrema.NbExt(), 1);
              double aMinimum = 1.0;
              for (int i = 1; i <= anExtrema.NbExt(); ++i)
              {
                aMinimum = std::min(aMinimum, std::sqrt(anExtrema.SquareDistance(i)));
                EXPECT_GE(anExtrema.Parameter(i), isTrimmed ? 0.1 * aSpan : 0.0);
                EXPECT_LE(anExtrema.Parameter(i), isTrimmed ? 0.8 * aSpan : 2 * aSpan);
              }
              EXPECT_NEAR(aMinimum, 1e-7, 4e-15);
            }
}
