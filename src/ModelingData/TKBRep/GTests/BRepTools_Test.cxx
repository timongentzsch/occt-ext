// Copyright (c) 2026 OPEN CASCADE SAS
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

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BinTools.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax3.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_Array1.hxx>
#include <OSD_Parallel.hxx>
#include <TopoDS_Shape.hxx>

#include <gtest/gtest.h>

#include <atomic>
#include <random>
#include <cmath>
#include <sstream>
#include <string>

namespace
{

// Write a box to a BRep string once, to be reused by all threads.
std::string MakeBoxBrepString()
{
  BRepPrimAPI_MakeBox aBox(10.0, 20.0, 30.0);
  EXPECT_FALSE(aBox.Shape().IsNull());

  std::ostringstream anOss;
  BRepTools::Write(aBox.Shape(), anOss);
  EXPECT_FALSE(anOss.str().empty());
  return anOss.str();
}

// Read a shape from a BRep buffer. Each call uses its own stream.
bool ReadFromBuffer(const std::string& theData, TopoDS_Shape& theShape)
{
  std::istringstream anIss(theData);
  BRep_Builder       aBuilder;
  BRepTools::Read(theShape, anIss, aBuilder);
  return !theShape.IsNull();
}

} // namespace

// Parallel BRep reading must produce shapes with identical volume properties.
// The underlying bug was that Strtod() shared a single Bigint memory pool
// across all threads, causing random mis-parses of numeric values in .brep
// files when read concurrently.
TEST(BRepTools_Test, ParallelRead_IdenticalMass)
{
  const std::string aBrepData = MakeBoxBrepString();

  // Single-threaded reference
  TopoDS_Shape aRefShape;
  ASSERT_TRUE(ReadFromBuffer(aBrepData, aRefShape));
  GProp_GProps aRefProps;
  BRepGProp::VolumeProperties(aRefShape, aRefProps);
  const double aRefMass = std::abs(aRefProps.Mass());
  ASSERT_GT(aRefMass, 0.0);

  constexpr int aNbReads = 200;

  NCollection_Array1<TopoDS_Shape> aShapes(0, aNbReads - 1);
  std::atomic<int>                 aFailCount{0};

  OSD_Parallel::For(0, aNbReads, [&](int theIndex) {
    TopoDS_Shape aShape;
    if (!ReadFromBuffer(aBrepData, aShape))
    {
      aFailCount.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    aShapes(theIndex) = aShape;
  });

  EXPECT_EQ(aFailCount.load(), 0) << "Some parallel reads failed";

  constexpr double aEps = 1.0e-6;
  for (int i = 0; i < aNbReads; ++i)
  {
    ASSERT_FALSE(aShapes(i).IsNull()) << "Shape " << i << " is null";
    GProp_GProps aProps;
    BRepGProp::VolumeProperties(aShapes(i), aProps);
    const double aMass = std::abs(aProps.Mass());
    EXPECT_NEAR(aMass, aRefMass, aEps)
      << "Shape " << i << " mass " << aMass << " != reference " << aRefMass;
  }
}

// Repeated parallel runs to catch intermittent failures.
TEST(BRepTools_Test, ParallelRead_RepeatedRuns)
{
  const std::string aBrepData = MakeBoxBrepString();

  TopoDS_Shape aRefShape;
  ASSERT_TRUE(ReadFromBuffer(aBrepData, aRefShape));
  GProp_GProps aRefProps;
  BRepGProp::VolumeProperties(aRefShape, aRefProps);
  const double aRefMass = std::abs(aRefProps.Mass());

  constexpr int    aNbReads = 100;
  constexpr int    aNbRuns  = 5;
  constexpr double aEps     = 1.0e-6;

  for (int aRun = 0; aRun < aNbRuns; ++aRun)
  {
    NCollection_Array1<TopoDS_Shape> aShapes(0, aNbReads - 1);
    std::atomic<int>                 aFailCount{0};

    OSD_Parallel::For(0, aNbReads, [&](int theIndex) {
      TopoDS_Shape aShape;
      if (!ReadFromBuffer(aBrepData, aShape))
      {
        aFailCount.fetch_add(1, std::memory_order_relaxed);
        return;
      }
      aShapes(theIndex) = aShape;
    });

    EXPECT_EQ(aFailCount.load(), 0) << "Run " << aRun << ": some reads failed";

    for (int i = 0; i < aNbReads; ++i)
    {
      ASSERT_FALSE(aShapes(i).IsNull()) << "Run " << aRun << ", shape " << i << " is null";
      GProp_GProps aProps;
      BRepGProp::VolumeProperties(aShapes(i), aProps);
      const double aMass = std::abs(aProps.Mass());
      EXPECT_NEAR(aMass, aRefMass, aEps) << "Run " << aRun << ", shape " << i << " mass mismatch";
    }
  }
}

// A shape read back from a file has the directions it was written with, bit for bit: the readers
// restore them with gp_Dir::Restored and gp_Dir2d::Restored. Lines, the normals
// of planes and lines on them in random oblique directions, through the binary and the text format.
// (A plane's X and Y directions are rebuilt from its normal and X direction by cross products when
// read, which this does not make exact.)
TEST(BRepTools_Test, ObliqueDirectionsReadBackBitForBit)
{
  std::mt19937_64                        aRandom(20261103u);
  std::uniform_real_distribution<double> aCoord(-1.0, 1.0);
  BRep_Builder                           aBuilder;
  TopoDS_Compound                        aCompound;
  aBuilder.MakeCompound(aCompound);
  for (int i = 0; i < 200; ++i)
  {
    const gp_Pnt anOrigin(10.0 * aCoord(aRandom), 10.0 * aCoord(aRandom), 10.0 * aCoord(aRandom));
    const gp_Dir aNormal(aCoord(aRandom), aCoord(aRandom), aCoord(aRandom));
    const gp_Dir anAlong = aNormal.Crossed(gp_Dir(aCoord(aRandom), aCoord(aRandom), aCoord(aRandom)));
    const occ::handle<Geom_Plane> aPlane = new Geom_Plane(gp_Ax3(anOrigin, aNormal, anAlong));
    TopoDS_Edge                   anEdge;
    aBuilder.MakeEdge(anEdge, new Geom_Line(anOrigin, anAlong), 1e-7);
    aBuilder.UpdateEdge(anEdge,
                        new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(aCoord(aRandom), aCoord(aRandom))),
                        aPlane,
                        TopLoc_Location(),
                        1e-7);
    aBuilder.Range(anEdge, 0.0, 1.0);
    aBuilder.Add(aCompound, anEdge);
  }
  std::ostringstream aBinary(std::ios::binary), aText;
  BinTools::Write(aCompound, aBinary);
  BRepTools::Write(aCompound, aText);
  TopoDS_Shape       aFromBinary, aFromText;
  std::istringstream aBinaryIn(aBinary.str(), std::ios::binary), aTextIn(aText.str());
  BinTools::Read(aFromBinary, aBinaryIn);
  BRepTools::Read(aFromText, aTextIn, aBuilder);
  for (const TopoDS_Shape& aRead : {aFromBinary, aFromText})
  {
    TopExp_Explorer aWritten(aCompound, TopAbs_EDGE), aBack(aRead, TopAbs_EDGE);
    int             aCount = 0;
    for (; aWritten.More() && aBack.More(); aWritten.Next(), aBack.Next(), ++aCount)
    {
      const TopoDS_Edge& aW = TopoDS::Edge(aWritten.Current());
      const TopoDS_Edge& aB = TopoDS::Edge(aBack.Current());
      double             aF, aL;
      const gp_Dir aDirW = occ::down_cast<Geom_Line>(BRep_Tool::Curve(aW, aF, aL))->Position().Direction();
      const gp_Dir aDirB = occ::down_cast<Geom_Line>(BRep_Tool::Curve(aB, aF, aL))->Position().Direction();
      ASSERT_EQ(aDirW.X(), aDirB.X());
      ASSERT_EQ(aDirW.Y(), aDirB.Y());
      ASSERT_EQ(aDirW.Z(), aDirB.Z());
      occ::handle<Geom2d_Curve> aPCurveW, aPCurveB;
      occ::handle<Geom_Surface> aSurfW, aSurfB;
      TopLoc_Location           aLocW, aLocB;
      BRep_Tool::CurveOnSurface(aW, aPCurveW, aSurfW, aLocW, aF, aL);
      BRep_Tool::CurveOnSurface(aB, aPCurveB, aSurfB, aLocB, aF, aL);
      const gp_Dir2d aD2W = occ::down_cast<Geom2d_Line>(aPCurveW)->Direction();
      const gp_Dir2d aD2B = occ::down_cast<Geom2d_Line>(aPCurveB)->Direction();
      ASSERT_EQ(aD2W.X(), aD2B.X());
      ASSERT_EQ(aD2W.Y(), aD2B.Y());
      const gp_Dir aNormalW = occ::down_cast<Geom_Plane>(aSurfW)->Position().Direction();
      const gp_Dir aNormalB = occ::down_cast<Geom_Plane>(aSurfB)->Position().Direction();
      ASSERT_EQ(aNormalW.X(), aNormalB.X());
      ASSERT_EQ(aNormalW.Y(), aNormalB.Y());
      ASSERT_EQ(aNormalW.Z(), aNormalB.Z());
    }
    EXPECT_EQ(aCount, 200);
    EXPECT_FALSE(aWritten.More() || aBack.More());
  }
}
