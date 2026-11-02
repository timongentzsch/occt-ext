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
// the product of the factors and the bounding box follows the factors exactly.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepPrimAPI_MakeCylinder.hxx>
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
