// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once

#include <vector>

namespace GeomAlgoExt
{

// Tensor Bernstein polynomials of degrees (n, m), coefficient (i, j) at i * (m + 1) + j: their
// products and degree raising, as hull bounds of surfaces and of their derivatives ask for.

// Binomial coefficients, from a table: products of Bernstein polynomials ask for them in their
// innermost loops.
inline double choose(int n, int k)
{
  static const auto table = [] {
    std::vector<std::vector<double>> rows(1, std::vector<double>{1});
    for (int i = 1; i <= 96; ++i)
    {
      std::vector<double> row(i + 1, 1);
      for (int j = 1; j < i; ++j)
        row[j] = rows[i - 1][j - 1] + rows[i - 1][j];
      rows.push_back(std::move(row));
    }
    return rows;
  }();
  if (n < int(table.size()))
    return table[n][k];
  double out = 1;
  for (int i = 1; i <= k; ++i)
    out = out * (n - k + i) / i;
  return out;
}

// The product of a, of degrees (n1, m1), and b, of degrees (n2, m2): of degrees (n1 + n2, m1 + m2).
template <class T>
std::vector<T> product(const std::vector<T>&      a,
                       int                        n1,
                       int                        m1,
                       const std::vector<double>& b,
                       int                        n2,
                       int                        m2)
{
  const int           n = n1 + n2, m = m1 + m2;
  std::vector<T>      out((n + 1) * (m + 1), T());
  std::vector<double> fv((m1 + 1) * (m2 + 1));
  for (int j = 0; j <= m1; ++j)
    for (int l = 0; l <= m2; ++l)
      fv[j * (m2 + 1) + l] = choose(m1, j) * choose(m2, l) / choose(m, j + l);
  for (int i = 0; i <= n1; ++i)
    for (int k = 0; k <= n2; ++k)
    {
      const double  cu  = choose(n1, i) * choose(n2, k) / choose(n, i + k);
      T*            row = &out[(i + k) * (m + 1)];
      const T*      x   = &a[i * (m1 + 1)];
      const double* y   = &b[k * (m2 + 1)];
      for (int j = 0; j <= m1; ++j)
        for (int l = 0; l <= m2; ++l)
          row[j + l] += x[j] * (y[l] * cu * fv[j * (m2 + 1) + l]);
    }
  return out;
}

// Raised to degrees (N, M): each coefficient a convex combination of the old ones.
template <class T>
std::vector<T> elevate(std::vector<T> a, int n, int m, int N, int M)
{
  for (; n < N; ++n)
  {
    std::vector<T> out((n + 2) * (m + 1), T());
    for (int i = 0; i <= n + 1; ++i)
      for (int j = 0; j <= m; ++j)
      {
        const double s     = double(i) / (n + 1);
        T            value = T();
        if (i > 0)
          value += a[(i - 1) * (m + 1) + j] * s;
        if (i <= n)
          value += a[i * (m + 1) + j] * (1 - s);
        out[i * (m + 1) + j] = value;
      }
    a = std::move(out);
  }
  for (; m < M; ++m)
  {
    std::vector<T> out((n + 1) * (m + 2), T());
    for (int i = 0; i <= n; ++i)
      for (int j = 0; j <= m + 1; ++j)
      {
        const double s     = double(j) / (m + 1);
        T            value = T();
        if (j > 0)
          value += a[i * (m + 1) + j - 1] * s;
        if (j <= m)
          value += a[i * (m + 1) + j] * (1 - s);
        out[i * (m + 2) + j] = value;
      }
    a = std::move(out);
  }
  return a;
}

} // namespace GeomAlgoExt
