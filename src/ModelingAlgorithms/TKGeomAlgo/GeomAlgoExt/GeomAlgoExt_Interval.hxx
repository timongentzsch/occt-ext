// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace GeomAlgoExt::interval
{
// Enclose exact operations on stored IEEE-754 doubles. Requires round-to-nearest,
// without fast-math/reassociation. Infinite intervals are unsupported: overflow
// declines the proof, and invalid intervals must never recover through clipping.
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53);

struct Range
{
  double low, high;
};

inline Range invalid()
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  return {nan, nan};
}

inline bool valid(Range x)
{
  return std::isfinite(x.low) && std::isfinite(x.high) && x.low <= x.high;
}

// Scalar rounding primitives do not certify their input; callers must check
// finite results before using them as bounds (as the interval operations do).
inline double roundDown(double x)
{
  return std::nextafter(x, -std::numeric_limits<double>::infinity());
}

inline double roundUp(double x)
{
  return std::nextafter(x, std::numeric_limits<double>::infinity());
}

inline Range point(double x)
{
  return std::isfinite(x) ? Range{x, x} : invalid();
}

namespace detail
{
inline Range enclose(double low, double high)
{
  if (!std::isfinite(low) || !std::isfinite(high) || low > high)
    return invalid();
  const Range result{roundDown(low), roundUp(high)};
  return valid(result) ? result : invalid();
}
} // namespace detail

inline Range add(Range a, Range b)
{
  if (!valid(a) || !valid(b))
    return invalid();
  return detail::enclose(a.low + b.low, a.high + b.high);
}

inline Range subtract(Range a, Range b)
{
  if (!valid(a) || !valid(b))
    return invalid();
  return detail::enclose(a.low - b.high, a.high - b.low);
}

inline Range multiply(Range a, Range b)
{
  if (!valid(a) || !valid(b))
    return invalid();
  const std::array<double, 4> products{a.low * b.low,
                                       a.low * b.high,
                                       a.high * b.low,
                                       a.high * b.high};
  for (double product : products)
    if (!std::isfinite(product))
      return invalid();
  return detail::enclose(*std::min_element(products.begin(), products.end()),
                         *std::max_element(products.begin(), products.end()));
}

inline Range dividePositive(Range a, Range b)
{
  if (!valid(a) || !valid(b) || !(b.low > 0))
    return invalid();
  const Range inverse = detail::enclose(1 / b.high, 1 / b.low);
  return multiply(a, inverse);
}

inline Range intersect(Range a, Range b)
{
  if (!valid(a) || !valid(b))
    return invalid();
  const Range result{std::max(a.low, b.low), std::min(a.high, b.high)};
  return valid(result) ? result : invalid();
}
} // namespace GeomAlgoExt::interval
