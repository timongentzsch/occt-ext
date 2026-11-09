// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once
#include <BSplCLib.hxx>
#include <NCollection_Array1.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Vec.hxx>
#include <vector>

namespace GeomAlgoExt
{
struct HermiteSite
{
  double parameter = 0;
  int    order = 0;
};

// The knot vector is expanded by multiplicity and clamped at both ends. Replace
// endpoint Greville sites by the requested jets (orders 0 through each maximum).
// Derivatives use the supplied knot parameter, without a length normalization.
inline std::vector<HermiteSite> hermiteSites(int                        degree,
                                             const std::vector<double>& flatKnots,
                                             int                        firstJetOrder,
                                             int                        lastJetOrder)
{
  const int count = int(flatKnots.size()) - degree - 1;
  if (degree < 1 || firstJetOrder < 0 || lastJetOrder < 0 || firstJetOrder > degree
      || lastJetOrder > degree || count < firstJetOrder + lastJetOrder + 2)
    return {};
  if (!std::is_sorted(flatKnots.begin(), flatKnots.end())
      || !std::all_of(flatKnots.begin(), flatKnots.end(), [](double x) { return std::isfinite(x); })
      || !(flatKnots.front() < flatKnots.back()))
    return {};
  for (int i = 0; i <= degree; ++i)
    if (flatKnots[i] != flatKnots.front()
        || flatKnots[flatKnots.size() - 1 - i] != flatKnots.back())
      return {};
  NCollection_Array1<double> knots(1, int(flatKnots.size())), sites(1, count);
  for (size_t i = 0; i < flatKnots.size(); ++i)
    knots(int(i) + 1) = flatKnots[i];
  BSplCLib::BuildSchoenbergPoints(degree, knots, sites);
  std::vector<HermiteSite> result;
  result.reserve(count);
  for (int i = 0; i < count; ++i)
  {
    if (i <= firstJetOrder)
      result.push_back({flatKnots.front(), i});
    else if (count - 1 - i <= lastJetOrder)
      result.push_back({flatKnots.back(), count - 1 - i});
    else
      result.push_back({sites(i + 1), 0});
  }
  return result;
}

// Values are the prescribed physical vectors or their parameter derivatives at
// the sites. On success they become the control vectors in the same spline
// basis. A singular collocation leaves the caller's values unchanged.
inline bool interpolateHermite(int                             degree,
                               const std::vector<double>&      flatKnots,
                               const std::vector<HermiteSite>& sites,
                               std::vector<gp_Vec>&            values)
{
  const int count = int(sites.size());
  if (count == 0 || values.size() != sites.size() || int(flatKnots.size()) != count + degree + 1)
    return false;
  NCollection_Array1<double> knots(1, int(flatKnots.size())), parameters(1, count);
  NCollection_Array1<int>    orders(1, count);
  std::vector<double>        poles(size_t(count) * 3);
  for (size_t i = 0; i < flatKnots.size(); ++i)
    knots(int(i) + 1) = flatKnots[i];
  for (int i = 0; i < count; ++i)
  {
    parameters(i + 1) = sites[i].parameter;
    orders(i + 1)     = sites[i].order;
    for (int j = 0; j < 3; ++j)
      poles[size_t(i) * 3 + j] = values[i].Coord(j + 1);
  }
  int singular = 0;
  BSplCLib::Interpolate(degree, knots, parameters, orders, 3, poles[0], singular);
  if (singular != 0)
    return false;
  for (int i = 0; i < count; ++i)
    values[i] = gp_Vec(poles[size_t(i) * 3], poles[size_t(i) * 3 + 1], poles[size_t(i) * 3 + 2]);
  return true;
}
} // namespace GeomAlgoExt
