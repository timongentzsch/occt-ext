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

// Chamfers and variable-radius fillets of BRepFilletExt on a straight box edge, where the removed
// material is a prism and has a closed form: its cross-section area times the edge length L.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepFilletExt_Blend.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_BezierCurve.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

// The edge of the box a x b x c that runs along z at (a, b).
TopoDS_Edge uprightAtFarCorner(const TopoDS_Shape& box, double a, double b)
{
  const auto at = [&](const double* l) {
    return near(l[0], a) && near(l[3], a) && near(l[1], b) && near(l[4], b);
  };
  const auto edges = edgesWhere(box, at);
  return edges.empty() ? TopoDS_Edge() : edges.front();
}

} // namespace

TEST(BRepFilletExt, LimitedBlendCarriesContactsAndSharpStops)
{
  for (bool circular : {false, true})
    for (bool variable : {false, true})
    {
      const TopoDS_Shape body = circular ? BRepPrimAPI_MakeCylinder(10, 20).Shape()
                                         : BRepPrimAPI_MakeBox(20, 20, 20).Shape();
      const auto candidates = circular
        ? edgesWhere(body, [](const double* b) { return b[2] > 19.999; })
        : std::vector<TopoDS_Edge>{uprightAtFarCorner(body, 20, 20)};
      ASSERT_FALSE(candidates.empty());
      const auto source = BRepExt::store(session(), body);
      const auto picked = BRepExt::store(session(), candidates.front());
      std::vector<BRepExt::blending::Join> joins;
      const auto output = BRepExt::blendRange(session(), source, picked, 1, 0, 0, .5, 0, 0,
        .2, .7, {}, BRepExt::failureHashBound,
        variable ? (circular ? std::vector<double>{0, 1, .5, 2, 1, 1}
                             : std::vector<double>{0, 1, 1, 2})
                 : std::vector<double>{}, &joins);
      const auto& made = BRepExt::get(session(), output.resultId);
      EXPECT_TRUE(BRepExt::invalidity(made).empty());
      NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> adjacent;
      TopExp::MapShapesAndAncestors(made, TopAbs_EDGE, TopAbs_FACE, adjacent);
      int contacts = 0, stops = 0;
      for (const auto& join : joins)
      {
        EXPECT_TRUE(adjacent.Contains(join.edge));
        if (join.role == BRepExt::blending::JoinRole::Contact) ++contacts;
        if (join.role == BRepExt::blending::JoinRole::Stop)
        { ++stops; EXPECT_EQ(join.order, 0); }
      }
      EXPECT_GE(contacts, 2);
      EXPECT_GE(stops, 2);
      for (int i = 1; i <= adjacent.Extent(); ++i)
      {
        bool constructed = false;
        for (const auto& face : adjacent.FindFromIndex(i))
        {
          bool old = false;
          for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
            old |= BRep_Tool::Surface(TopoDS::Face(face)) == BRep_Tool::Surface(TopoDS::Face(f.Current()));
          constructed |= !old;
        }
        if (!constructed) continue;
        bool recorded = false;
        for (const auto& join : joins) recorded |= join.edge.IsSame(adjacent.FindKey(i));
        EXPECT_TRUE(recorded) << "unclassified constructed edge " << i;
      }
    }
}

// An equal-distance chamfer of distance d removes the triangular prism of area d^2 / 2: a wrong
// setback on either face, or a bevel that is not planar, changes the volume by far more than the
// tolerance.
TEST(BRepFilletExt, EqualChamferRemovesATriangularPrism)
{
  const double       a = 6, b = 4, c = 3;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  for (const double d : {0.3, 1.0})
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
    const auto     result =
      BRepExt::chamfer(session(), body, {picked}, d, d, 0, {}, BRepExt::failureHashBound, 0, false);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c - 0.5 * d * d * c, 1e-9 * a * b * c) << "d=" << d;
  }
}

// An asymmetric chamfer takes d1 along the reference face and d2 along the other: the prism
// removed has area d1 d2 / 2, whichever of the two faces is the reference.
TEST(BRepFilletExt, AsymmetricChamferRemovesTheRightTriangle)
{
  const double       a = 6, b = 4, c = 3, d1 = 1.5, d2 = 0.5;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    faces;
  TopExp::MapShapesAndAncestors(box, TopAbs_EDGE, TopAbs_FACE, faces);
  const auto& adjacent = faces.FindFromKey(edge);
  ASSERT_EQ(adjacent.Extent(), 2);
  for (const TopoDS_Shape& reference : adjacent)
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge),
                   ref         = BRepExt::store(session(), reference);
    const auto          result = BRepExt::chamfer(session(),
                                         body,
                                                  {picked},
                                         d1,
                                         d2,
                                         ref,
                                                  {},
                                         BRepExt::failureHashBound,
                                         0,
                                         false);
    const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c - 0.5 * d1 * d2 * c, 1e-9 * a * b * c);
  }
}

// A variable fillet whose law is constant is the circular fillet: the quarter-disc complement of
// area r^2 (1 - pi/4). It exercises the variable-radius path (law points, sections along the
// spine) on a case whose answer does not depend on it.
TEST(BRepFilletExt, ConstantLawVariableFilletIsTheCircularFillet)
{
  const double       a = 6, b = 4, c = 3, r = 0.5;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  const auto     result    = BRepExt::filletVariable(session(),
                                              body,
                                                     {picked},
                                                     {3},
                                                     {0, r, 0.5, r, 1, r},
                                                     {},
                                              BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c - r * r * (1 - M_PI / 4) * c, 1e-6 * r * r * c);
}

// A law that grows along the edge removes more material than its smaller end alone and less than
// its larger end alone, since every cross-section of the blend lies between the two circular
// ones: the volume lies between the two constant-radius closed forms. A law ignored, or a radius
// taken from the wrong parameter range, leaves it outside.
TEST(BRepFilletExt, GrowingLawRemovesBetweenTheTwoConstantRadiusVolumes)
{
  const double       a = 6, b = 4, c = 3, r0 = 0.3, r1 = 0.9;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  const auto     result    = BRepExt::filletVariable(session(),
                                              body,
                                                     {picked},
                                                     {2},
                                                     {0, r0, 1, r1},
                                                     {},
                                              BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const double removed = a * b * c - volume(made), k = (1 - M_PI / 4) * c;
  EXPECT_GT(removed, k * r0 * r0);
  EXPECT_LT(removed, k * r1 * r1);
}

// A chamfer, like every operation, refuses what it cannot do with a typed failure instead of
// returning a wrong body: an asymmetric chamfer without a reference face has no side to take d1 on.
TEST(BRepFilletExt, AsymmetricChamferWithoutReferenceIsRefused)
{
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(6, 4, 3).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, 6, 4);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  EXPECT_THROW(BRepExt::chamfer(session(),
                                body,
                                {picked},
                                1.5,
                                0.5,
                                0,
                                {},
                                BRepExt::failureHashBound,
                                0,
                                false),
               std::exception);
}

// A chamfer, and one whose setback grows linearly along the edge, kept between the planes square
// to the edge at 0.2 and 0.7 of it: the prism removed there only, of section d^2 / 2, or of
// d(z)^2 / 2 integrated over the range for d = d0 + (d1 - d0) z / c. Ranges for every section a
// blend takes go through the blend of the edge alone (blendRange).
TEST(BRepFilletExt, ChamferInARangeRemovesThePrismThereOnly)
{
  const double       a = 6, b = 4, c = 10, d0 = 0.5, d1 = 1.5, z0 = 0.2 * c, z1 = 0.7 * c;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  // The integral of (d0 + k z)^2 / 2 over [z0, z1], k the law's slope.
  const double k        = (d1 - d0) / c;
  const auto   integral = [&](double z) { return std::pow(d0 + k * z, 3) / (6 * k); };
  for (const bool grows : {false, true})
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
    const auto     result =
      BRepExt::blendRange(session(),
                          body,
                          picked,
                          d0,
                          0,
                          2,
                          0.5,
                          0,
                          0,
                          0.2,
                          0.7,
                          {},
                          BRepExt::failureHashBound,
                          grows ? std::vector<double>{0, d0, 1, d1} : std::vector<double>{});
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    const double removed = grows ? integral(z1) - integral(z0) : 0.5 * d0 * d0 * (z1 - z0);
    EXPECT_NEAR(volume(made), a * b * c - removed, 1e-9 * a * b * c) << (grows ? "grows" : "");
  }
}

// A blend kept between the planes square to the edge at 0.17 and 0.63 of it, on a box edge and on
// a cylinder's rim, convex and in a cavity, conic and chamfer (G2 takes the section policy of
// BRepBlendExt_Policy), at two scales. The
// sections are the same all along, so the change within the range is that share of the change the
// blend of the whole edge makes (of length along the box edge, of angle round the rim; volumes of
// NURBS copies, as OpenCascade's volume of extrusion faces is off by 1e-3 here), and no vertex of
// the result is looser than the body's.
TEST(BRepFilletExt, LimitedBlendTakesTheRangeShareOfTheChange)
{
  const auto exactVolume = [](const TopoDS_Shape& shape) {
    GProp_GProps properties;
    BRepGProp::VolumeProperties(BRepBuilderAPI_NurbsConvert(shape).Shape(), properties, 1e-12);
    return properties.Mass();
  };
  for (bool rim : {false, true})
    for (bool cavity : {false, true})
      for (int law : {0, 2})
        for (double scale : {1., 100.})
        {
          SCOPED_TRACE(testing::Message() << "rim " << rim << " cavity " << cavity << " law " << law
                                          << " scale " << scale);
          TopoDS_Shape body = rim ? BRepPrimAPI_MakeCylinder(10, 20).Shape()
                                  : BRepPrimAPI_MakeBox(20, 20, 20).Shape();
          const auto onEdge = [&](const double* b) {
            return rim ? b[2] > 19.999 && b[5] < 20.001
                       : near(b[0], 0) && near(b[3], 0) && near(b[1], 0) && near(b[4], 0) && b[5] - b[2] > 19;
          };
          if (cavity)
            body = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(gp_Pnt(-30, -30, -10), gp_Pnt(30, 30, 30)).Shape(), body)
                     .Shape();
          const auto edges = edgesWhere(body, onEdge);
          ASSERT_EQ(edges.size(), 1u);
          gp_Trsf scaling;
          scaling.SetScale(gp_Pnt(0, 0, 0), scale);
          BRepBuilderAPI_Transform transform(body, scaling, true);
          body             = transform.Shape();
          const auto edge  = transform.ModifiedShape(edges[0]);
          double loosest = 0;
          for (TopExp_Explorer v(body, TopAbs_VERTEX); v.More(); v.Next())
            loosest = std::max(loosest, BRep_Tool::Tolerance(TopoDS::Vertex(v.Current())));
          const uint32_t id = BRepExt::store(session(), body), picked = BRepExt::store(session(), edge);
          const double   before = exactVolume(body);
          const auto     whole  = BRepExt::blendRange(session(), id, picked, scale, 0, law, .5, 0, 0, 0, 1, {},
                                                 BRepExt::failureHashBound);
          const auto     ranged = BRepExt::blendRange(session(), id, picked, scale, 0, law, .5, 0, 0, .17, .63, {},
                                                  BRepExt::failureHashBound);
          const TopoDS_Shape& made = BRepExt::get(session(), ranged.resultId);
          EXPECT_TRUE(BRepExt::invalidity(made).empty());
          const double change = exactVolume(BRepExt::get(session(), whole.resultId)) - before;
          EXPECT_NEAR(exactVolume(made) - before, (.63 - .17) * change, 1e-6 * std::abs(change));
          for (TopExp_Explorer v(made, TopAbs_VERTEX); v.More(); v.Next())
            EXPECT_LE(BRep_Tool::Tolerance(TopoDS::Vertex(v.Current())), loosest);
        }
}

// A chamfer at 30 degrees off a reference face whose setback on it grows along the edge, over the
// whole edge and over a range: on a right angle the other setback is d tan 30, and the section
// removed d^2 tan 30 / 2, integrated as above.
TEST(BRepFilletExt, AngleChamferWhoseSetbackGrowsRemovesItsSections)
{
  const double       a = 6, b = 4, c = 10, d0 = 0.5, d1 = 1.5, angle = M_PI / 6;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  const auto faces =
    facesWhere(box, [&](const double* l) { return near(l[0], a) && near(l[3], a); });
  ASSERT_EQ(faces.size(), 1u);
  const double k        = (d1 - d0) / c;
  const auto   integral = [&](double z) {
    return std::tan(angle) * std::pow(d0 + k * z, 3) / (6 * k);
  };
  for (const auto& [start, end] : {std::pair{0.0, 1.0}, std::pair{0.2, 0.7}})
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge),
                   reference   = BRepExt::store(session(), faces[0]);
    const auto          result = BRepExt::blendRange(session(),
                                            body,
                                            picked,
                                            d0,
                                            0,
                                            2,
                                            0.5,
                                            reference,
                                            angle,
                                            start,
                                            end,
                                                     {},
                                            BRepExt::failureHashBound,
                                                     {0, d0, 1, d1});
    const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made),
                a * b * c - (integral(end * c) - integral(start * c)),
                1e-9 * a * b * c)
      << start << " " << end;
  }
}

// An open parabolic rim has both ends in y>5 and its middle outside it.
// The retained portions are two runs, although there are only two crossings.
TEST(BRepFilletExt, OpenLimitedContourCannotRetainTwoRuns)
{
  for (double level : {1.0, 5.0, 9.0})
  {
  SCOPED_TRACE(level);
  NCollection_Array1<gp_Pnt> poles(1, 3);
  poles(1) = gp_Pnt(0, 10, 0); poles(2) = gp_Pnt(10, -10, 0); poles(3) = gp_Pnt(20, 10, 0);
  const Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(poles);
  BRepBuilderAPI_MakeWire wire;
  wire.Add(BRepBuilderAPI_MakeEdge(curve).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(20, 10, 0), gp_Pnt(20, 20, 0)).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(20, 20, 0), gp_Pnt(0, 20, 0)).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0, 20, 0), gp_Pnt(0, 10, 0)).Edge());
  const TopoDS_Shape body = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wire.Wire()).Face(), gp_Vec(0,0,10));
  const auto selected = edgesWhere(body, [](const double* b) { return near(b[2],10) && near(b[5],10) && b[1]<1; });
  ASSERT_EQ(selected.size(), 1);
  const auto id = BRepExt::store(session(), body), edge = BRepExt::store(session(), selected[0]);
  try
  {
    BRepExt::filletStations(session(), id, {edge}, .5, {}, {}, {edge},
      {10 * (1 - std::sqrt(level / 10)), level, 10, 0, 1, 0}, {}, BRepExt::failureHashBound);
    ADD_FAILURE() << "an open limited contour retained two disconnected runs";
  }
  catch (const BRepExt::Failure& failure)
  {
    EXPECT_EQ(failure.reason, BRepExt::unsupported_geometry) << failure.what();
  }
}
}
