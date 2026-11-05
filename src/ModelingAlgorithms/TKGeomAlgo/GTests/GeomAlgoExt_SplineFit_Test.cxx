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

#include <GeomAlgoExt_PatchFit.hxx>
#include <gtest/gtest.h>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_PatchCompact.hxx>
#include <Message_ProgressIndicator.hxx>
#include <GeomAlgoExt_Progress.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <Precision.hxx>
#include <cmath>

// A wrong basis offset or derivative row fails this closed-form Bernstein oracle.
TEST(GeomAlgoExtSplineFit, CubicBasisJets)
{
  const std::vector<double> knots{0, 0, 0, 0, 1, 1, 1, 1};
  for (double u : {0.0, .125, .5, .875, 1.0})
  {
    const auto   b = BRepExt::splineBasis("basis", knots, 3, u);
    const double v = 1 - u;
    const double value[4]{v * v * v, 3 * u * v * v, 3 * u * u * v, u * u * u};
    const double first[4]{-3 * v * v, 3 * v * v - 6 * u * v, 6 * u * v - 3 * u * u, 3 * u * u};
    const double second[4]{6 * v, 18 * u - 12, 6 - 18 * u, 6 * u};
    const double third[4]{-6, 18, -18, 6};
    EXPECT_EQ(b.first, 0);
    for (int i = 0; i < 4; ++i)
    {
      EXPECT_NEAR(b.value[0][i], value[i], Precision::PConfusion());
      EXPECT_NEAR(b.value[1][i], first[i], Precision::PConfusion());
      EXPECT_NEAR(b.value[2][i], second[i], Precision::PConfusion());
      EXPECT_NEAR(b.value[3][i], third[i], Precision::PConfusion());
    }
  }
}

// Tridiagonal SPD matrices: the expected x is independent of the factor/solve.
TEST(GeomAlgoExtSplineFit, BandedCholeskyResidual)
{
  for (bool paired : {false, true})
    for (int n : {4, 17, 65})
    {
      GeomAlgoExt::Band   a(n, 1, paired);
      std::vector<double> x(n), expected(n);
      for (int i = 0; i < n; ++i)
        expected[i] = std::sin(double(i));
      for (int i = 0; i < n; ++i)
      {
        a.add(i, i, 2);
        if (i)
          a.add(i, i - 1, -.25);
        x[i] = 2 * expected[i] - (i ? .25 * expected[i - 1] : 0)
               - (i + 1 < n ? .25 * expected[i + 1] : 0);
      }
      ASSERT_TRUE(a.factor());
      a.solve(x);
      for (int i = 0; i < n; ++i)
        EXPECT_NEAR(x[i], expected[i], Precision::PConfusion());
    }
}

// A cubic tensor polynomial is uniquely determined by its 4x4 values.
// Crossed U/V strides, incorrect fixed-pole elimination or non-finite results fail.
TEST(GeomAlgoExtSplineFit, TensorPolynomialInterpolation)
{
  for (double scale : {.01, 1., 100.})
  {
    BRepExt::SplineNet net;
    net.degreeU = net.degreeV = 3;
    net.knotsU = net.knotsV = {0, 0, 0, 0, 1, 1, 1, 1};
    net.poles.resize(16, gp_XYZ(0, 0, 0));
    auto exact = [&](double u, double v) {
      return gp_Pnt(scale * (u + v * v),
                    scale * (v - u * u),
                    scale * (u * u * u - 2 * u * v + v * v * v));
    };
    std::vector<BRepExt::FitCondition> conditions;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
      {
        const double u = i / 3., v = j / 3.;
        const auto   p = exact(u, v);
        for (int c = 0; c < 3; ++c)
        {
          BRepExt::FitCondition row{};
          row.u        = u;
          row.v        = v;
          row.count    = 1;
          row.terms[0] = {0, 0, gp_XYZ(c == 0, c == 1, c == 2)};
          row.value    = p.Coord(c + 1);
          row.weight   = 1;
          row.within   = Precision::Confusion() * scale;
          conditions.push_back(row);
        }
      }
    BRepExt::fitNet("interpolation", net, std::vector<char>(16, 0), conditions, {0, 0}, {});
    auto surface = BRepExt::surfaceOf(net);
    for (int i = 0; i <= 16; ++i)
      for (int j = 0; j <= 16; ++j)
        EXPECT_LE(surface->Value(i / 16., j / 16.).Distance(exact(i / 16., j / 16.)),
                  Precision::Confusion() * scale);
  }
}

namespace
{
BRepExt::SplineNet planeNet(double height)
{
  BRepExt::SplineNet net;
  net.degreeU = net.degreeV = 3;
  net.knotsU = net.knotsV = {0, 0, 0, 0, 1, 1, 1, 1};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      net.poles.emplace_back(i / 3., j / 3., height);
  return net;
}

class StopIndicator : public Message_ProgressIndicator
{
  bool UserBreak() override { return true; }

  void Show(const Message_ProgressScope&, const bool) override {}
};
} // namespace

// The squared second/third derivative energy is nonnegative. An affine plane
// has zero energy and the four fixed corners determine it uniquely. This checks
// both coordinate layouts and directional movement, with nonuniform knot spans
// and fixed-pole elimination, against that independent minimum.
TEST(GeomAlgoExtSplineFit, FairFitRecoversAnAffinePlaneFromPerturbedPoles)
{
  for (int degree : {3, 5})
    for (double scale : {.01, 1., 100.})
      for (int mode = 0; mode < 3; ++mode)
      {
        SCOPED_TRACE(testing::Message() << degree << ", " << scale << ", " << mode);
        BRepExt::SplineNet net;
        net.degreeU = net.degreeV = degree;
        net.knotsU.assign(degree + 1, 0);
        net.knotsU.insert(net.knotsU.end(), {.21, .63});
        net.knotsU.insert(net.knotsU.end(), degree + 1, 1);
        net.knotsV.assign(degree + 1, 0);
        net.knotsV.insert(net.knotsV.end(), {.37});
        net.knotsV.insert(net.knotsV.end(), degree + 1, 1);
        const int nu = net.countU(), nv = net.countV();
        const auto greville = [degree](const std::vector<double>& knots, int i) {
          double sum = 0;
          for (int k = 1; k <= degree; ++k)
            sum += knots[i + k];
          return sum / degree;
        };
        const auto exact = [scale](double u, double v) {
          return gp_XYZ(scale * (1 + 2 * u - v),
                        scale * (-2 + u + 3 * v),
                        scale * (3 - u + 2 * v));
        };
        const gp_XYZ direction = gp_XYZ(1, 2, 3).Normalized();
        std::vector<char> fixed(nu * nv, 0);
        for (int i = 0; i < nu; ++i)
          for (int j = 0; j < nv; ++j)
          {
            const bool corner = (i == 0 || i == nu - 1) && (j == 0 || j == nv - 1);
            gp_XYZ p = exact(greville(net.knotsU, i), greville(net.knotsV, j));
            if (!corner)
              p += direction * (scale * std::sin(1.7 * i + .8 * j));
            net.poles.push_back(p);
            fixed[i * nv + j] = corner;
          }
        std::vector<BRepExt::FitCondition> conditions;
        if (mode == 1)
        {
          BRepExt::FitCondition row{};
          row.u = .31;
          row.v = .59;
          row.count = 1;
          row.terms[0] = {0, 0, direction};
          row.value = exact(row.u, row.v).Dot(direction);
          row.weight = 1;
          row.within = Precision::Confusion() * scale;
          conditions.push_back(row);
        }
        BRepExt::fitNet("affine", net, fixed, conditions, {1, .1}, {},
                        mode == 2 ? &direction : nullptr);
        const auto surface = BRepExt::surfaceOf(net);
        for (int i = 0; i <= 8; ++i)
          for (int j = 0; j <= 8; ++j)
            EXPECT_LE(surface->Value(i / 8., j / 8.).Distance(gp_Pnt(exact(i / 8., j / 8.))),
                      Precision::Confusion() * scale);
      }
}

// A wrong normal-system RHS or pole update fails the defining plane distances.
TEST(GeomAlgoExtSplineFit, RimCorrectionToParallelPlane)
{
  auto       net    = planeNet(0);
  const auto places = [](double) {
    std::vector<BRepExt::RimPlace> out;
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; j <= 8; ++j)
        out.push_back({gp_Pnt(i / 8., j / 8., .125), gp_XY(i / 8., j / 8.), gp_XY(0, 0), false});
    return out;
  };
  const double error =
    BRepExt::correctRim("rim", net, places, gp_XYZ(0, 0, 1), Precision::Confusion(), 0, 16);
  EXPECT_LE(error, Precision::Confusion());
  auto surface = BRepExt::surfaceOf(net);
  for (int i = 0; i <= 16; ++i)
    for (int j = 0; j <= 16; ++j)
      EXPECT_NEAR(surface->Value(i / 16., j / 16.).Z(), .125, Precision::Confusion());
}

// A false interval certificate would accept displacement greater than its bound.
TEST(GeomAlgoExtSplineFit, KnotRemovalBoundOnAffineSurface)
{
  auto surface = BRepExt::surfaceOf(planeNet(.125));
  surface->InsertUKnot(.5, 1, Precision::PConfusion());
  surface->InsertVKnot(.5, 1, Precision::PConfusion());
  auto reduced = BRepExt::compactPatchSurface("compact", surface, Precision::Confusion());
  ASSERT_TRUE(reduced.has_value());
  for (int i = 0; i <= 32; ++i)
    for (int j = 0; j <= 32; ++j)
    {
      const auto point = reduced->surface->Value(i / 32., j / 32.);
      EXPECT_LE(point.Distance(gp_Pnt(i / 32., j / 32., .125)), Precision::Confusion());
      EXPECT_LE(point.Distance(surface->Value(i / 32., j / 32.)), reduced->displacement);
    }
}

// Optional optimization catches OCCT builder failures, but must never swallow a stop.
TEST(GeomAlgoExtSplineFit, CompactCancellationEscapesOptionalFailureCatch)
{
  auto surface = BRepExt::surfaceOf(planeNet(0));
  surface->InsertUKnot(.5, 1, Precision::PConfusion());
  Handle(Message_ProgressIndicator) stop  = new StopIndicator;
  const auto                        range = stop->Start();
  BRepExt::ProgressScope            progress(range);
  try
  {
    BRepExt::compactPatchSurface("compact", surface, Precision::Confusion());
    FAIL() << "cancel did not escape";
  }
  catch (const BRepExt::Failure& failure)
  {
    EXPECT_EQ(failure.reason, "cancelled");
    EXPECT_EQ(failure.message, "Operation cancelled");
  }
  for (int i = 0; i <= 8; ++i)
    for (int j = 0; j <= 8; ++j)
      EXPECT_LE(surface->Value(i / 8., j / 8.).Distance(gp_Pnt(i / 8., j / 8., 0)),
                Precision::Confusion());
}

// A leaked or unrestored thread-local progress range either cancels the inner
// certificate or lets the outer cancelled call complete. Exercise real proposals.
TEST(GeomAlgoExtSplineFit, NestedNativeProgressRestoresCancellation)
{
  for (double height : {-.125, 0., .25})
  {
    const auto make = [](double z) {
      auto result = BRepExt::surfaceOf(planeNet(z));
      result->InsertUKnot(.5, 1, Precision::PConfusion());
      return result;
    };
    Handle(Message_ProgressIndicator) stop = new StopIndicator;
    {
      BRepExt::ProgressScope outer(stop->Start());
      {
        BRepExt::ProgressScope inner(Message_ProgressRange{});
        auto                   input = make(height);
        auto reduced = BRepExt::compactPatchSurface("nested", input, Precision::Confusion());
        ASSERT_TRUE(reduced.has_value());
        for (int i = 0; i <= 16; ++i)
          for (int j = 0; j <= 16; ++j)
            EXPECT_LE(
              reduced->surface->Value(i / 16., j / 16.).Distance(input->Value(i / 16., j / 16.)),
              reduced->displacement);
      }
      EXPECT_THROW(BRepExt::compactPatchSurface("outer", make(height), Precision::Confusion()),
                   BRepExt::Failure);
    }
    EXPECT_NO_THROW(BRepExt::compactPatchSurface("restored", make(height), Precision::Confusion()));
  }
}
