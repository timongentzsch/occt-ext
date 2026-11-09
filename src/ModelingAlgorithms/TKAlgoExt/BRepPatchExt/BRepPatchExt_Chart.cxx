// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepPatchExt_Chart.hxx>
#include <GeomAlgoExt_Bernstein.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <Geom_BezierSurface.hxx>

#include <algorithm>
#include <array>
#include <cmath>

namespace BRepExt
{
namespace
{
const char* const op = "patch";

//! The point of the triangle (a, b, c) nearest to p, by its weights on the corners (Ericson,
//! Real-Time Collision Detection, 5.1.5: the Voronoi regions of the corners, the edges and the
//! inside in turn).
std::array<double, 3> nearestWeights(const gp_Pnt& p,
                                     const gp_Pnt& a,
                                     const gp_Pnt& b,
                                     const gp_Pnt& c)
{
  const gp_Vec ab(a, b), ac(a, c), ap(a, p);
  const double d1 = ab.Dot(ap), d2 = ac.Dot(ap);
  if (d1 <= 0 && d2 <= 0)
    return {1, 0, 0};
  const gp_Vec bp(b, p);
  const double d3 = ab.Dot(bp), d4 = ac.Dot(bp);
  if (d3 >= 0 && d4 <= d3)
    return {0, 1, 0};
  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0)
  {
    const double v = d1 / (d1 - d3);
    return {1 - v, v, 0};
  }
  const gp_Vec cp(c, p);
  const double d5 = ab.Dot(cp), d6 = ac.Dot(cp);
  if (d6 >= 0 && d5 <= d6)
    return {0, 0, 1};
  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0)
  {
    const double w = d2 / (d2 - d6);
    return {1 - w, 0, w};
  }
  const double va = d3 * d6 - d5 * d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
  {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return {0, 1 - w, w};
  }
  const double denominator = 1 / (va + vb + vc);
  const double v = vb * denominator, w = vc * denominator;
  return {1 - v - w, v, w};
}

//! A polynomial of degrees (n, m) in Bernstein form, one vector per coefficient, cut at the
//! middle of u (`alongU`) or of v by de Casteljau: the two halves, each over the whole square.
std::array<std::vector<gp_XYZ>, 2> halves(const std::vector<gp_XYZ>& c, int n, int m, bool alongU)
{
  std::array<std::vector<gp_XYZ>, 2> out{c, c};
  const int                          lines = alongU ? m + 1 : n + 1, count = alongU ? n + 1 : m + 1;
  for (int line = 0; line < lines; ++line)
  {
    const auto at = [&](int k) {
      return alongU ? static_cast<size_t>(k) * (m + 1) + line
                    : static_cast<size_t>(line) * (m + 1) + k;
    };
    std::vector<gp_XYZ> work(count);
    for (int k = 0; k < count; ++k)
      work[k] = c[at(k)];
    for (int r = 1; r < count; ++r)
    {
      out[0][at(r - 1)]     = work[0];
      out[1][at(count - r)] = work[count - r];
      for (int k = 0; k + r < count; ++k)
        work[k] = (work[k] + work[k + 1]) * 0.5;
    }
    out[0][at(count - 1)] = work[0];
    out[1][at(0)]         = work[0];
  }
  return out;
}

//! The least cosine of the coefficients against their mean direction, cut in four up to `depth`
//! times where it is not positive; nothing where it is not positive at the last.
std::optional<double> coneMargin(const std::vector<gp_XYZ>& c, int n, int m, int depth)
{
  gp_XYZ mean(0, 0, 0);
  for (const gp_XYZ& v : c)
    mean += v;
  double least = -1;
  if (mean.Modulus() > 0)
  {
    mean /= mean.Modulus();
    least = 1;
    for (const gp_XYZ& v : c)
    {
      const double length = v.Modulus();
      least               = std::min(least, length > 0 ? v.Dot(mean) / length : -1.0);
    }
  }
  if (least > 0)
    return least;
  if (depth == 0)
    return std::nullopt;
  double worst = 1;
  for (const auto& half : halves(c, n, m, true))
    for (const auto& quarter : halves(half, n, m, false))
    {
      const auto margin = coneMargin(quarter, n, m, depth - 1);
      if (!margin)
        return std::nullopt;
      worst = std::min(worst, *margin);
    }
  return worst;
}

//! How often a cell whose coefficients are not on one side is cut in four before it counts as
//! not shown regular: four times, cells of a 16th of a span each way. Cutting converges to the
//! polynomial's own values (the control net of a halved Bernstein polynomial approaches it as
//! 4^-depth), so a cell that is regular is shown so at some depth; this is how far it is tried.
constexpr int coneCuts = 4;

//! How often the fair disc is laid out again as an even lattice in its own chart and faired
//! there (see proxyChart). A judgement, not derived: measured by how little the chart moves on
//! the next (see the probe in dox/extensions/operations/BRepPatchExt.md).
constexpr int relays = 2;

//! The nodes that triangles of `nodes` keep, of which the first `frontier` are the boundary in
//! turn, renumbered in their order: the chart's `source` and `triangles`.
ConformalChart meshed(const std::vector<gp_XY>& nodes, size_t frontier)
{
  ConformalChart chart;
  chart.triangles = triangulateInside(nodes, frontier);
  std::vector<int> kept(nodes.size(), -1);
  for (const auto& t : chart.triangles)
    for (int k : t)
      kept[k] = 0;
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (kept[i] < 0)
      continue;
    kept[i] = static_cast<int>(chart.source.size());
    chart.source.push_back(nodes[i]);
  }
  for (size_t k = 0; k < frontier; ++k)
    if (kept[k] != static_cast<int>(k))
      fail(op, no_solution, "the loop's disc leaves out a place of the loop");
  for (auto& t : chart.triangles)
    for (int& k : t)
      k = kept[k];
  return chart;
}

//! An even triangular lattice inside `boundary`, of the boundary's mean step, kept half a step
//! off it; the boundary first.
ConformalChart lattice(const std::vector<gp_XY>& boundary)
{
  const size_t m      = boundary.size();
  double       length = 0;
  gp_XY        low = boundary.front(), high = low;
  for (size_t k = 0; k < m; ++k)
  {
    length += (boundary[(k + 1) % m] - boundary[k]).Modulus();
    low.SetCoord(std::min(low.X(), boundary[k].X()), std::min(low.Y(), boundary[k].Y()));
    high.SetCoord(std::max(high.X(), boundary[k].X()), std::max(high.Y(), boundary[k].Y()));
  }
  const double       step = length / m, rowStep = step * std::sqrt(3.0) / 2;
  std::vector<gp_XY> nodes = boundary;
  for (int row = 0; low.Y() + row * rowStep < high.Y(); ++row)
    for (double x = low.X() + ((row & 1) ? step / 2 : 0); x < high.X(); x += step)
    {
      const gp_XY p(x, low.Y() + row * rowStep);
      bool        in      = false;
      double      nearest = HUGE_VAL;
      for (size_t i = 0, j = m - 1; i < m; j = i++)
      {
        const gp_XY edge = boundary[i] - boundary[j], to = p - boundary[j];
        nearest = std::min(
          nearest,
          (to - edge * std::clamp(to.Dot(edge) / edge.SquareModulus(), 0.0, 1.0)).Modulus());
        if ((boundary[i].Y() > p.Y()) != (boundary[j].Y() > p.Y())
            && p.X() < boundary[j].X() + edge.X() * to.Y() / edge.Y())
          in = !in;
      }
      if (in && nearest > step / 2)
        nodes.push_back(p);
    }
  return meshed(nodes, m);
}

//! A symmetric positive definite matrix by rows from the first entry that is not zero to the
//! diagonal (a skyline), and its Cholesky factor in place. The fairing's own, of internal linkage:
//! an inline solver shared with the fits would let the linker keep this unit's copy of it for
//! theirs as well.
struct Skyline
{
  std::vector<int>                 first;
  std::vector<std::vector<double>> rows;

  explicit Skyline(int size)
      : first(size),
        rows(size)
  {
    for (int i = 0; i < size; ++i)
    {
      first[i] = i;
      rows[i].assign(1, 0);
    }
  }

  double& at(int r, int c) { return rows[r][c - first[r]]; }

  //! Adds to entry (r, c), r >= c.
  void add(int r, int c, double value)
  {
    if (c < first[r])
    {
      rows[r].insert(rows[r].begin(), first[r] - c, 0);
      first[r] = c;
    }
    at(r, c) += value;
  }

  bool factor()
  {
    for (int i = 0; i < static_cast<int>(rows.size()); ++i)
      for (int j = first[i]; j <= i; ++j)
      {
        double sum = at(i, j);
        for (int k = std::max(first[i], first[j]); k < j; ++k)
          sum -= at(i, k) * at(j, k);
        if (j < i)
          at(i, j) = sum / at(j, j);
        else if (!(sum > 0) || !std::isfinite(sum))
          return false;
        else
          at(i, i) = std::sqrt(sum);
      }
    return true;
  }

  void solve(std::vector<double>& x)
  {
    const int n = static_cast<int>(rows.size());
    for (int i = 0; i < n; ++i)
    {
      for (int k = first[i]; k < i; ++k)
        x[i] -= at(i, k) * x[k];
      x[i] /= at(i, i);
    }
    for (int i = n - 1; i >= 0; --i)
    {
      x[i] /= at(i, i);
      for (int k = first[i]; k < i; ++k)
        x[k] -= at(i, k) * x[i];
    }
  }
};

//! A place of a triangle, by its weights on the corners, that is to be at `point`.
struct Pin
{
  std::array<int, 3>    nodes;
  std::array<double, 3> weights;
  gp_Pnt                point;
};

//! What a pin counts for against the umbrella of one node, both lengths. A judgement, not
//! derived: enough that the disc passes within a small share of the guides' height of them
//! (see the probe in dox/extensions/operations/BRepPatchExt.md), not so much that it creases at them.
constexpr double pinWeight = 100;

//! The triangle of `chart` that holds `uv` in its source, or the nearest, and the weights of
//! its corners there.
std::pair<std::array<int, 3>, std::array<double, 3>> locate(const ConformalChart& chart,
                                                            const gp_XY&          uv)
{
  double                best = -HUGE_VAL;
  std::array<int, 3>    nodes{};
  std::array<double, 3> weights{};
  for (const auto& t : chart.triangles)
  {
    const gp_XY &a = chart.source[t[0]], &b = chart.source[t[1]], &c = chart.source[t[2]];
    const double whole = (b - a).Crossed(c - a);
    if (!(std::abs(whole) > 0))
      continue;
    const std::array<double, 3> w{(b - uv).Crossed(c - uv) / whole,
                                  (c - uv).Crossed(a - uv) / whole,
                                  (a - uv).Crossed(b - uv) / whole};
    const double                least = std::min({w[0], w[1], w[2]});
    if (least > best)
    {
      best    = least;
      nodes   = t;
      weights = w;
    }
  }
  // Just outside, the nearest triangle's weights are clamped onto it.
  double sum = 0;
  for (double& w : weights)
    sum += (w = std::max(w, 0.0));
  for (double& w : weights)
    w /= sum;
  return {nodes, weights};
}

//! The places in space of `chart`'s nodes, its first nodes being `rim`'s places (see proxyChart).
void fair(ConformalChart& chart, const std::vector<ChartRim>& rim, const std::vector<Pin>& pins)
{
  const int m = static_cast<int>(rim.size()), n = static_cast<int>(chart.source.size());
  // Neighbours by the triangles, and for the loop's nodes the outer ring's two beside them:
  // outer node k stands beyond the piece from loop node k to the next.
  std::vector<std::vector<int>> neighbours(n);
  for (const auto& t : chart.triangles)
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (a != b)
          neighbours[t[a]].push_back(t[b]);
  for (auto& list : neighbours)
  {
    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
  }
  for (int k = 0; k < m; ++k)
  {
    neighbours[k].push_back(n + (k + m - 1) % m);
    neighbours[k].push_back(n + k);
  }
  // Every node's place, known or not; and its number among the unknowns, -1 where known.
  std::vector<gp_XYZ> place(n + m, gp_XYZ(0, 0, 0));
  std::vector<int>    unknown(n + m, -1);
  std::vector<gp_XY>  disc(n + m);
  for (int i = 0; i < n; ++i)
    disc[i] = chart.source[i];
  // An outer node is numbered as if it stood at the middle of its piece.
  for (int k = 0; k < m; ++k)
  {
    place[k]            = rim[k].point.XYZ();
    const gp_Pnt& next  = rim[(k + 1) % m].point;
    disc[n + k]         = (chart.source[k] + chart.source[(k + 1) % m]) / 2;
    const double height = rim[k].point.Distance(next) * std::sqrt(3.0) / 2;
    if (rim[k].leave.SquareMagnitude() > 0)
      place[n + k] =
        (rim[k].point.XYZ() + next.XYZ()) / 2 - rim[k].leave.Normalized().XYZ() * height;
  }
  // Numbered by where they lie in the disc, row after row, so that the normal equations have
  // a narrow band: each umbrella reaches a node's neighbours only.
  std::vector<int> order;
  for (int i = m; i < n; ++i)
    order.push_back(i);
  for (int k = 0; k < m; ++k)
    if (rim[k].leave.SquareMagnitude() == 0)
      order.push_back(n + k);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return disc[a].Y() != disc[b].Y() ? disc[a].Y() < disc[b].Y() : disc[a].X() < disc[b].X();
  });
  for (size_t i = 0; i < order.size(); ++i)
    unknown[order[i]] = static_cast<int>(i);
  const int count = static_cast<int>(order.size());
  if (count == 0)
    fail(op, no_solution, "the loop's disc has nothing inside");

  // The normal equations of the umbrellas: per node, its own place less the mean of its
  // neighbours'; and of the pins, each a place of a triangle that is to be at its point.
  struct Term
  {
    int    at;
    double weight;
  };

  struct Row
  {
    std::vector<Term> terms;
    gp_XYZ            known = gp_XYZ(0, 0, 0);
    double            weight;
  };

  std::vector<Row> rows;
  const auto       row = [&](const std::vector<std::pair<int, double>>& at, double weight) {
    Row out{{}, gp_XYZ(0, 0, 0), weight};
    for (const auto& [j, w] : at)
      if (unknown[j] >= 0)
        out.terms.push_back({unknown[j], w});
      else
        out.known += place[j] * w;
    rows.push_back(std::move(out));
  };
  for (int i = 0; i < n; ++i)
  {
    std::vector<std::pair<int, double>> at{{i, -1}};
    for (int j : neighbours[i])
      at.push_back({j, 1.0 / neighbours[i].size()});
    row(at, 1);
  }
  for (const Pin& pin : pins)
  {
    row({{pin.nodes[0], pin.weights[0]},
         {pin.nodes[1], pin.weights[1]},
         {pin.nodes[2], pin.weights[2]}},
        pinWeight);
    rows.back().known -= pin.point.XYZ();
  }
  Skyline                            normal(count);
  std::array<std::vector<double>, 3> right;
  for (auto& r : right)
    r.assign(count, 0);
  for (const Row& r : rows)
    for (const Term& a : r.terms)
    {
      for (const Term& b : r.terms)
        if (b.at <= a.at)
          normal.add(a.at, b.at, r.weight * a.weight * b.weight);
      for (int c = 0; c < 3; ++c)
        right[c][a.at] -= r.weight * a.weight * r.known.Coord(c + 1);
    }
  if (!normal.factor())
    fail(op, no_solution, "the loop's disc has no fairest places");
  for (auto& r : right)
    normal.solve(r);
  chart.point.resize(n);
  for (int i = 0; i < n; ++i)
    chart.point[i] = unknown[i] < 0
                       ? gp_Pnt(place[i])
                       : gp_Pnt(right[0][unknown[i]], right[1][unknown[i]], right[2][unknown[i]]);
}

} // namespace

// The places in space of the disc are the fairest that a lattice can have where they span the
// loop: those that make the sum of the squares of the umbrella Laplacian, at every node of the
// loop and inside it, least (discrete thin-plate fairing: Kobbelt, "Discrete fairing", 1997;
// Botsch and Kobbelt, "An intuitive framework for real-time freeform modeling", 2004). The
// loop's own places are held. A fill that is to be tangent to a face leaves its edge as the
// face would go on across it, and is held to that by a ring of nodes outside the disc, one per
// piece of the loop, standing in the face: at the piece's middle less `leave` times the height
// of an even triangle on it. With the umbrella of a node of the loop then near nothing, the
// nodes inside next to it stand where the face's tangent plane goes on across the loop (the
// discrete form of the condition that fixes the first derivative across the boundary, as two
// rings of fixed nodes do for the biharmonic equation in Botsch and Kobbelt). Where nothing is
// asked, the outer node is free, and the loop's umbrella is then as the fit makes it: a natural
// end.
// The umbrella is the Laplacian of the lattice's own parameters, which are not the surface's:
// first the unit disc, with the loop at even steps of its length on its rim; then the chart
// that the fair disc lays flat to, in which an even lattice is near even on the surface as
// well, so that its umbrella is near the surface's own Laplacian (the cotangent Laplacian of a
// mesh of even triangles: Pinkall and Polthier 1993). Faired and laid flat again on that
// lattice, `relays` times. From the first of those on, the disc is pinned to the points the
// fill is to pass through (the guides'), each where the disc before was nearest it.
ConformalChart proxyChart(const std::vector<ChartRim>& rim, const std::vector<gp_Pnt>& through)
{
  PatchPhase phase("proxy-chart");
  const int  m = static_cast<int>(rim.size());
  if (m < 3)
    fail(op, invalid_input, "the loop to be charted has fewer than three places");
  std::vector<double> along(m + 1, 0);
  for (int k = 0; k < m; ++k)
    along[k + 1] = along[k] + rim[k].point.Distance(rim[(k + 1) % m].point);
  const double length = along[m];
  if (!(length > 0))
    fail(op, invalid_input, "the loop to be charted has no length");
  std::vector<gp_XY> circle;
  for (int k = 0; k < m; ++k)
  {
    const double angle = 2 * M_PI * along[k] / length;
    circle.emplace_back(std::cos(angle), std::sin(angle));
  }
  ConformalChart chart = lattice(circle);
  fair(chart, rim, {});
  layFlat(chart);
  for (int relay = 0; relay < relays; ++relay)
  {
    ConformalChart next = lattice(std::vector<gp_XY>(chart.flat.begin(), chart.flat.begin() + m));
    // The points the fill is to pass through, where the disc nearest them lies in the new
    // lattice's source, the chart before.
    std::vector<Pin> pins;
    for (const gp_Pnt& point : through)
    {
      double     distance = 0;
      const auto at       = locate(next, onChart(chart, point, distance));
      pins.push_back({at.first, at.second, point});
    }
    fair(next, rim, pins);
    layFlat(next);
    chart = std::move(next);
  }
  return chart;
}

std::vector<double> nodeBends(const ConformalChart& chart)
{
  const size_t        n = chart.point.size();
  std::vector<gp_XYZ> laplace(n, gp_XYZ(0, 0, 0));
  std::vector<double> area(n, 0), angle(n, 0);
  std::vector<char>   onBoundary(n, 0);
  for (const int b : chart.boundary)
    onBoundary[b] = 1;
  for (const auto& t : chart.triangles)
  {
    const double whole = gp_Vec(chart.point[t[0]], chart.point[t[1]])
                           .Crossed(gp_Vec(chart.point[t[0]], chart.point[t[2]]))
                           .Magnitude()
                         / 2;
    for (int k = 0; k < 3; ++k)
    {
      const int    i = t[k], j = t[(k + 1) % 3], l = t[(k + 2) % 3];
      const gp_Vec a(chart.point[l], chart.point[i]), b(chart.point[l], chart.point[j]);
      const double cross = a.Crossed(b).Magnitude();
      if (!(cross > 0))
        continue;
      const double cotangent = a.Dot(b) / cross;
      angle[l] += std::atan2(cross, a.Dot(b));
      area[l] += whole / 3;
      const gp_XYZ edge = chart.point[i].XYZ() - chart.point[j].XYZ();
      laplace[i] += edge * (cotangent / 2);
      laplace[j] -= edge * (cotangent / 2);
    }
  }
  // No bend tighter than the triangles round a node is read from them: one over their mean edge.
  std::vector<double> edges(n, 0);
  std::vector<int>    counted(n, 0);
  for (const auto& t : chart.triangles)
    for (int k = 0; k < 3; ++k)
    {
      const double length = chart.point[t[k]].Distance(chart.point[t[(k + 1) % 3]]);
      for (const int end : {t[k], t[(k + 1) % 3]})
      {
        edges[end] += length;
        ++counted[end];
      }
    }
  std::vector<double> bend(n, 0);
  for (size_t i = 0; i < n; ++i)
  {
    if (onBoundary[i] || !(area[i] > 0))
      continue;
    const double mean  = laplace[i].Modulus() / (2 * area[i]),
                 gauss = (2 * M_PI - angle[i]) / area[i];
    bend[i] = std::min(mean + std::sqrt(std::max(0.0, mean * mean - gauss)), counted[i] / edges[i]);
  }
  for (const auto& t : chart.triangles)
    for (int k = 0; k < 3; ++k)
      for (int e = 1; e < 3; ++e)
        if (onBoundary[t[k]] && !onBoundary[t[(k + e) % 3]])
          bend[t[k]] = std::max(bend[t[k]], bend[t[(k + e) % 3]]);
  return bend;
}

std::vector<double> nodeStretch(const ConformalChart& chart)
{
  const size_t        n = chart.point.size();
  std::vector<double> flat(n, 0), space(n, 0);
  for (const auto& t : chart.triangles)
  {
    const double inSpace = gp_Vec(chart.point[t[0]], chart.point[t[1]])
                             .Crossed(gp_Vec(chart.point[t[0]], chart.point[t[2]]))
                             .Magnitude();
    const double laid =
      std::abs((chart.flat[t[1]] - chart.flat[t[0]]).Crossed(chart.flat[t[2]] - chart.flat[t[0]]));
    for (const int k : t)
    {
      flat[k] += laid;
      space[k] += inSpace;
    }
  }
  std::vector<double> out(n, 1);
  for (size_t i = 0; i < n; ++i)
    if (space[i] > 0)
      out[i] = std::sqrt(flat[i] / space[i]);
  return out;
}

gp_XY onChart(const ConformalChart& chart, const gp_Pnt& point, double& distance)
{
  distance = HUGE_VAL;
  gp_XY best(0, 0);
  for (const auto& t : chart.triangles)
  {
    const auto   w = nearestWeights(point, chart.point[t[0]], chart.point[t[1]], chart.point[t[2]]);
    const gp_Pnt foot(chart.point[t[0]].XYZ() * w[0] + chart.point[t[1]].XYZ() * w[1]
                      + chart.point[t[2]].XYZ() * w[2]);
    const double d = foot.Distance(point);
    if (d < distance)
    {
      distance = d;
      best     = chart.flat[t[0]] * w[0] + chart.flat[t[1]] * w[1] + chart.flat[t[2]] * w[2];
    }
  }
  return best;
}

// The search is over the directions of one half of the sphere, laid out by the golden angle
// (a Fibonacci lattice: Gonzalez, "Measurement of areas on a sphere using Fibonacci and
// latitude-longitude lattices", Math. Geosci. 42, 2010), with the mean plane's own normal
// among them; the best of them is then made exact for the signs it gives the normals (below).
// Every answer is a direction's own least, a bound from below of the best.
double graphMargin(const std::vector<gp_Vec>& normals, const gp_Vec (&axes)[3])
{
  if (normals.empty())
    return 1;
  const auto least = [&](const gp_XYZ& d) {
    const gp_Vec direction = axes[0] * d.X() + axes[1] * d.Y() + axes[2] * d.Z();
    double       out       = 1;
    for (const gp_Vec& n : normals)
      out = std::min(out, std::abs(n.Dot(direction)) / n.Magnitude());
    return out;
  };
  constexpr int directions = 4096;
  gp_XYZ        best(0, 0, 1);
  double        margin = least(best);
  for (int k = 0; k < directions; ++k)
  {
    const double z = 1 - (k + 0.5) / directions, r = std::sqrt(1 - z * z),
                 phi = k * M_PI * (3 - std::sqrt(5.0));
    const gp_XYZ d(r * std::cos(phi), r * std::sin(phi), z);
    const double here = least(d);
    if (here > margin)
    {
      margin = here;
      best   = d;
    }
  }
  // Near the best direction every normal keeps its sign along it, and with the signs fixed the
  // best is exact: the largest over unit d of the least s n . d is the distance from the origin to
  // the convex hull of the signed normals s n, where it does not hold the origin, d pointing at its
  // nearest point (the duality of the support function: Boyd and Vandenberghe, Convex
  // Optimization, 2004, 8.1.1). That point is found by Gilbert's iteration (IEEE J. Robotics and
  // Automation 2, 1966: step to the nearest point of the segment to the vertex furthest back),
  // which leaves a bound from below at every step.
  std::vector<gp_XYZ> faced;
  const gp_Vec        towards = axes[0] * best.X() + axes[1] * best.Y() + axes[2] * best.Z();
  for (const gp_Vec& n : normals)
    faced.push_back((n.Dot(towards) < 0 ? -n : n).Normalized().XYZ());
  gp_XYZ nearest = faced.front();
  for (const gp_XYZ& a : faced)
    if (a.Dot(towards.XYZ()) < nearest.Dot(towards.XYZ()))
      nearest = a;
  for (int step = 0; step < 4096; ++step)
  {
    const gp_XYZ* back = &faced.front();
    for (const gp_XYZ& a : faced)
      if (a.Dot(nearest) < back->Dot(nearest))
        back = &a;
    const double length = nearest.SquareModulus();
    if (!(length > 0) || length - back->Dot(nearest) <= 1e-15 * length)
      break;
    const gp_XYZ way = *back - nearest;
    nearest += way * std::clamp(-nearest.Dot(way) / way.SquareModulus(), 0.0, 1.0);
  }
  if (nearest.Modulus() > 0)
  {
    // The least over the signed normals along the direction found, which is a margin whatever
    // the iteration has reached.
    double here = 1;
    for (const gp_XYZ& a : faced)
      here = std::min(here, a.Dot(nearest) / nearest.Modulus());
    margin = std::max(margin, here);
  }
  return margin;
}

// A Bezier patch of degrees (p, q) over a cell has S_u = p / du sum (P[i+1][j] - P[i][j])
// B^(p-1)_i B^q_j and S_v likewise, so each component of S_u x S_v is a sum of products of
// Bernstein polynomials of degrees (p - 1, q) and (p, q - 1), one of degrees (2p - 1, 2q - 1)
// (Farin, Curves and Surfaces for CAGD, 5th ed., 2002, sections 5.8 and 17). By the convex
// hull property its values on the cell are convex combinations of its coefficients: where
// these all have a positive part along one direction, so has every value, and none is zero.
// This is the surface form of the Bernstein test of a Jacobian's sign used for planar
// parametrisations (Gravesen, Evgrafov, Nguyen, Nortoft, "Planar parametrization in
// isogeometric analysis", 2014).
std::optional<double> regularMargin(const occ::handle<Geom_BSplineSurface>& surface,
                                    const std::function<bool(double, double, double, double)>& used)
{
  PatchPhase                                phase("regularity");
  GeomConvert_BSplineSurfaceToBezierSurface cut(surface);
  const int                                 cellsU = cut.NbUPatches(), cellsV = cut.NbVPatches();
  NCollection_Array2<occ::handle<Geom_BezierSurface>> patches(1, cellsU, 1, cellsV);
  cut.Patches(patches);
  NCollection_Array1<double> knotsU(1, cellsU + 1), knotsV(1, cellsV + 1);
  cut.UKnots(knotsU);
  cut.VKnots(knotsV);
  double worst = 1;
  for (int a = 1; a <= cellsU; ++a)
    for (int b = 1; b <= cellsV; ++b)
    {
      if (!used(knotsU(a), knotsU(a + 1), knotsV(b), knotsV(b + 1)))
        continue;
      const auto& patch = patches(a, b);
      const int   p = patch->UDegree(), q = patch->VDegree();
      // The derivatives' coefficients, per component: S_u of degrees (p - 1, q), S_v of (p, q - 1).
      std::array<std::vector<double>, 3> su, sv;
      for (int c = 0; c < 3; ++c)
      {
        su[c].assign(static_cast<size_t>(p) * (q + 1), 0);
        sv[c].assign(static_cast<size_t>(p + 1) * q, 0);
      }
      for (int i = 0; i <= p; ++i)
        for (int j = 0; j <= q; ++j)
          for (int c = 0; c < 3; ++c)
          {
            const double here = patch->Pole(i + 1, j + 1).Coord(c + 1);
            if (i < p)
              su[c][i * (q + 1) + j] = patch->Pole(i + 2, j + 1).Coord(c + 1) - here;
            if (j < q)
              sv[c][i * q + j] = patch->Pole(i + 1, j + 2).Coord(c + 1) - here;
          }
      // The scale of the cell changes no sign and no cosine.
      std::vector<gp_XYZ> normal(static_cast<size_t>(2 * p) * (2 * q), gp_XYZ(0, 0, 0));
      for (int c = 0; c < 3; ++c)
      {
        const int  d = (c + 1) % 3, e = (c + 2) % 3;
        const auto plus  = GeomAlgoExt::product(su[d], p - 1, q, sv[e], p, q - 1);
        const auto minus = GeomAlgoExt::product(su[e], p - 1, q, sv[d], p, q - 1);
        for (size_t k = 0; k < normal.size(); ++k)
          normal[k].SetCoord(c + 1, plus[k] - minus[k]);
      }
      const auto margin = coneMargin(normal, 2 * p - 1, 2 * q - 1, coneCuts);
      if (!margin)
        return std::nullopt;
      worst = std::min(worst, *margin);
    }
  return worst;
}

} // namespace BRepExt
