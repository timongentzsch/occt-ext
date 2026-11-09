// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// CSLib_Class2d keeps large polygons in a bucket index over their V range. The index may only skip
// edges that cannot change the answer, so every query must give what a scan of the whole polygon
// gives. The reference below is that scan, written out from the definition of the classification:
// a point is uncertain when a vertex lies within the tolerance box or an edge passes within the V
// tolerance; otherwise it is inside when a ray to +U crosses an odd number of edges, and uncertain
// again when a corner of the tolerance box falls on the other side.

#include <CSLib_Class2d.hxx>

#include <NCollection_Array1.hxx>
#include <Precision.hxx>
#include <gp_Pnt2d.hxx>

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

namespace
{

using Pnt = std::pair<double, double>;

bool crossings(const std::vector<Pnt>& thePolygon, const double thePx, const double thePy)
{
  int          aCount = 0;
  const size_t aN     = thePolygon.size();
  for (size_t i = 0; i < aN; ++i)
  {
    const Pnt&   aA  = thePolygon[i];
    const Pnt&   aB  = thePolygon[(i + 1) % aN];
    const double aAx = aA.first - thePx, aAy = aA.second - thePy;
    const double aBx = aB.first - thePx, aBy = aB.second - thePy;
    if ((aAy < 0.0) == (aBy < 0.0))
    {
      continue;
    }
    if (aAx > 0.0 && aBx > 0.0)
    {
      ++aCount;
    }
    else if (aAx > 0.0 || aBx > 0.0)
    {
      if (aAx - aAy * (aBx - aAx) / (aBy - aAy) > 0.0)
      {
        ++aCount;
      }
    }
  }
  return (aCount & 1) != 0;
}

CSLib_Class2d::Result reference(const std::vector<Pnt>& thePolygon,
                                const double            theTolU,
                                const double            theTolV,
                                const double            thePx,
                                const double            thePy)
{
  const size_t aN = thePolygon.size();
  for (size_t i = 0; i < aN; ++i)
  {
    const Pnt&   aA  = thePolygon[i];
    const Pnt&   aB  = thePolygon[(i + 1) % aN];
    const double aBx = aB.first - thePx, aBy = aB.second - thePy;
    if (aBx < theTolU && aBx > -theTolU && aBy < theTolV && aBy > -theTolV)
    {
      return CSLib_Class2d::Result_Uncertain;
    }
    const double aEdgeDx = aB.first - aA.first;
    if ((aA.first - thePx) * aBx < 0.0 && std::abs(aEdgeDx) > Precision::PConfusion())
    {
      const double aY = aB.second - (aB.second - aA.second) / aEdgeDx * aBx;
      if (aY - thePy >= -theTolV && aY - thePy <= theTolV)
      {
        return CSLib_Class2d::Result_Uncertain;
      }
    }
  }
  const bool isInside = crossings(thePolygon, thePx, thePy);
  if (theTolU > 0.0 || theTolV > 0.0)
  {
    for (const double aSx : {-theTolU, theTolU})
    {
      for (const double aSy : {-theTolV, theTolV})
      {
        if (isInside != crossings(thePolygon, thePx + aSx, thePy + aSy))
        {
          return CSLib_Class2d::Result_Uncertain;
        }
      }
    }
  }
  return isInside ? CSLib_Class2d::Result_Inside : CSLib_Class2d::Result_Outside;
}

// A closed curve r(t) around (0.5, 0.5) with n vertices, rough enough to have many edges per
// height and several crossings of a horizontal line.
std::vector<Pnt> wavyLoop(std::mt19937& theRandom,
                          const int     theN,
                          const double  theWaves,
                          const double  theAmplitude)
{
  std::uniform_real_distribution<double> aNoise(-0.002, 0.002);
  std::vector<Pnt>                       aPolygon;
  for (int i = 0; i < theN; ++i)
  {
    const double aT = 2.0 * M_PI * i / theN;
    const double aR = 0.3 + theAmplitude * std::sin(theWaves * aT) + aNoise(theRandom);
    aPolygon.emplace_back(0.5 + aR * std::cos(aT), 0.5 + aR * std::sin(aT));
  }
  return aPolygon;
}

CSLib_Class2d build(const std::vector<Pnt>& thePolygon, const double theTolU, const double theTolV)
{
  NCollection_Array1<gp_Pnt2d> aPoints(1, static_cast<int>(thePolygon.size()));
  for (size_t i = 0; i < thePolygon.size(); ++i)
  {
    aPoints.SetValue(static_cast<int>(i) + 1, gp_Pnt2d(thePolygon[i].first, thePolygon[i].second));
  }
  return CSLib_Class2d(aPoints, theTolU, theTolV, 0.0, 0.0, 1.0, 1.0);
}

} // namespace

// The index is built for polygons of 32 vertices and more; the answers must be those of the scan
// at every size on both sides of that, with random points, points on the vertices, points a
// fraction of the tolerance off a vertex, and points on the middle of an edge, so that the
// uncertain answers (which come from the edges nearest the point) are compared as well as the
// inside and outside ones.
TEST(CSLib_Class2dTest, IndexedAnswersEqualTheScanOfTheWholePolygon)
{
  std::mt19937                           aRandom(20260102);
  std::uniform_real_distribution<double> aUnit(0.0, 1.0);
  int                                    aUncertain = 0, aInside = 0, aOutside = 0;
  for (const int aN : {12, 31, 32, 33, 64, 257, 1000, 4096})
  {
    for (const double aTol : {0.0, 1e-7, 1e-4, 3e-3})
    {
      const std::vector<Pnt> aPolygon = wavyLoop(aRandom, aN, 7.0, 0.08);
      const CSLib_Class2d    aClass   = build(aPolygon, aTol, aTol);
      std::vector<Pnt>       aQueries;
      for (int k = 0; k < 400; ++k)
      {
        aQueries.emplace_back(aUnit(aRandom), aUnit(aRandom));
      }
      for (size_t i = 0; i < aPolygon.size(); i += 1 + aPolygon.size() / 60)
      {
        const Pnt& aA = aPolygon[i];
        const Pnt& aB = aPolygon[(i + 1) % aPolygon.size()];
        aQueries.push_back(aA);
        aQueries.emplace_back(aA.first + 0.4 * aTol, aA.second - 0.6 * aTol);
        aQueries.emplace_back(aA.first + 1.5 * aTol, aA.second + 1.5 * aTol);
        aQueries.emplace_back((aA.first + aB.first) / 2, (aA.second + aB.second) / 2);
        aQueries.emplace_back((aA.first + aB.first) / 2 + 0.5 * aTol, (aA.second + aB.second) / 2);
      }
      for (const Pnt& aQ : aQueries)
      {
        const CSLib_Class2d::Result aGot = aClass.SiDans(gp_Pnt2d(aQ.first, aQ.second));
        CSLib_Class2d::Result       aWant;
        if (aQ.first < -aTol || aQ.first > 1.0 + aTol || aQ.second < -aTol
            || aQ.second > 1.0 + aTol)
        {
          aWant = CSLib_Class2d::Result_Outside;
        }
        else
        {
          aWant = reference(aPolygon, aTol, aTol, aQ.first, aQ.second);
        }
        ASSERT_EQ(aGot, aWant) << "n=" << aN << " tol=" << aTol << " at " << aQ.first << ","
                               << aQ.second;
        aUncertain += aGot == CSLib_Class2d::Result_Uncertain;
        aInside += aGot == CSLib_Class2d::Result_Inside;
        aOutside += aGot == CSLib_Class2d::Result_Outside;
      }
    }
  }
  // The comparison covers all three answers, not only the common one.
  EXPECT_GT(aUncertain, 100);
  EXPECT_GT(aInside, 100);
  EXPECT_GT(aOutside, 100);
}

// Horizontal edges, repeated heights and a polygon pressed against the bounds are where a bucket
// by height goes wrong: a staircase has every edge at one of a few heights, and its vertices lie
// exactly on bucket boundaries.
TEST(CSLib_Class2dTest, IndexedAnswersEqualTheScanOnAStaircaseOnBucketBoundaries)
{
  std::vector<Pnt> aPolygon;
  const int        aSteps = 40;
  for (int i = 0; i < aSteps; ++i)
  {
    aPolygon.emplace_back(double(i) / aSteps, double(i) / aSteps);
    aPolygon.emplace_back(double(i + 1) / aSteps, double(i) / aSteps);
  }
  aPolygon.emplace_back(1.0, 1.0);
  aPolygon.emplace_back(0.0, 1.0);
  const double        aTol   = 1e-6;
  const CSLib_Class2d aClass = build(aPolygon, aTol, aTol);
  for (int i = 0; i <= 4 * aSteps; ++i)
  {
    for (int j = 0; j <= 4 * aSteps; ++j)
    {
      const double aX = double(i) / (4 * aSteps), aY = double(j) / (4 * aSteps);
      for (const double aShift : {0.0, 1e-7, -1e-7, 2e-6})
      {
        ASSERT_EQ(aClass.SiDans(gp_Pnt2d(aX + aShift, aY)),
                  reference(aPolygon, aTol, aTol, aX + aShift, aY))
          << aX << "," << aY << " shift " << aShift;
      }
    }
  }
}
