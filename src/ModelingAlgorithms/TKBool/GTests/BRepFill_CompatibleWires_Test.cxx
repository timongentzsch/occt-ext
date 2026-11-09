// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepFill_CompatibleWires.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Ellipse.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <TopLoc_Location.hxx>
#include <Precision.hxx>
#include <gp_Trsf.hxx>
#include <cmath>
#include <gtest/gtest.h>

// Every new ellipse cut must be a radial hit of a polygon vertex. Poses and
// scales exercise plane locations and the conversion to plane coordinates.
TEST(BRepFill_CompatibleWires, PlanarPolarCutsMeetEllipses)
{
  for (int seed = 0; seed < 9; ++seed)
  {
    const int n = 3 + seed % 5;
    const double scale = std::pow(10.0, seed % 3 - 1);
    BRepBuilderAPI_MakeWire polygon;
    for (int i = 0; i < n; ++i)
    {
      const double t0 = 2*M_PI*i/n, t1 = 2*M_PI*(i+1)/n;
      polygon.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(3*scale*std::cos(t0), 3*scale*std::sin(t0), 0),
                                         gp_Pnt(3*scale*std::cos(t1), 3*scale*std::sin(t1), 0)));
    }
    const auto ellipse = new Geom_Ellipse(gp_Ax2(gp_Pnt(0,0,4*scale), gp_Dir(0,0,1)),
                                           4*scale, 2*scale);
    TopoDS_Wire ring = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(ellipse)).Wire();
    gp_Trsf pose;
    pose.SetRotation(gp_Ax1(gp_Pnt(2,-3,1), gp_Dir(0.3,-0.5,1)), 0.31*seed);
    TopLoc_Location location(pose);
    NCollection_Sequence<TopoDS_Shape> sections;
    sections.Append(polygon.Wire().Moved(location));
    sections.Append(ring.Moved(location));
    BRepFill_CompatibleWires compatible(sections);
    compatible.Perform(false);
    ASSERT_TRUE(compatible.IsDone()) << seed;
    int counts[2] = {0,0};
    for (int row = 0; row < 2; ++row)
      for (TopExp_Explorer e(compatible.Shape().Value(row+1), TopAbs_EDGE); e.More(); e.Next())
        ++counts[row];
    EXPECT_EQ(counts[0], counts[1]);
    EXPECT_GE(counts[1], n);
    for (TopExp_Explorer v(compatible.Shape().Value(2), TopAbs_VERTEX); v.More(); v.Next())
    {
      gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(v.Current()));
      p.Transform(pose.Inverted());
      EXPECT_NEAR(p.Z(), 4*scale, Precision::Confusion());
      EXPECT_NEAR(std::pow(p.X()/(4*scale),2) + std::pow(p.Y()/(2*scale),2), 1, 1e-9);
    }
  }
}
