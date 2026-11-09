// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepSurfaceExt_Polyspline.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>

#include <BRepAlgoAPI_Check.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>

namespace BRepExt
{
namespace
{
constexpr const char* op = "polySpline";
using Quad               = std::array<int, 4>;
using Key                = std::pair<int, int>;

Key key(int a, int b)
{
  return std::minmax(a, b);
}

struct Mesh
{
  std::vector<gp_XYZ> p;
  std::vector<Quad>   q;
  std::vector<char>   corner;
};

struct Edge
{
  int              a, b;
  std::vector<int> faces;
};

struct Topology
{
  std::vector<Edge>               e;
  std::map<Key, int>              ids;
  std::vector<std::vector<int>>   faces, edges, boundary;
  std::vector<std::array<int, 4>> sides;

  explicit Topology(const Mesh& m)
      : faces(m.p.size()),
        edges(m.p.size()),
        boundary(m.p.size())
  {
    for (size_t f = 0; f < m.q.size(); ++f)
    {
      std::array<int, 4> sidesOf{};
      std::set<int>      unique(m.q[f].begin(), m.q[f].end());
      if (unique.size() != 4)
        fail(op, invalid_input, "a quad repeats a vertex");
      for (int i = 0; i < 4; ++i)
      {
        int a = m.q[f][i], b = m.q[f][(i + 1) % 4];
        faces[a].push_back(static_cast<int>(f));
        auto inserted = ids.emplace(key(a, b), static_cast<int>(e.size()));
        if (inserted.second)
          e.push_back({a, b, {}});
        auto& edge = e[inserted.first->second];
        if (!edge.faces.empty() && edge.a == a)
          fail(op, invalid_input, "quad winding is inconsistent across an edge");
        edge.faces.push_back(static_cast<int>(f));
        if (edge.faces.size() > 2)
          fail(op, invalid_input, "nonmanifold edge (more than two faces)");
        sidesOf[i] = inserted.first->second;
      }
      sides.push_back(sidesOf);
    }
    for (size_t i = 0; i < e.size(); ++i)
    {
      auto& edge = e[i];
      edges[edge.a].push_back(static_cast<int>(i));
      edges[edge.b].push_back(static_cast<int>(i));
      if (edge.faces.size() == 1)
      {
        boundary[edge.a].push_back(edge.b);
        boundary[edge.b].push_back(edge.a);
      }
    }
    for (size_t v = 0; v < m.p.size(); ++v)
    {
      if (faces[v].empty())
        fail(op, invalid_input, "unused cage vertex");
      if (!boundary[v].empty() && boundary[v].size() != 2)
        fail(op, invalid_input, "nonmanifold boundary vertex");
      std::set<int>   visited;
      std::queue<int> todo;
      todo.push(faces[v][0]);
      while (!todo.empty())
      {
        int f = todo.front();
        todo.pop();
        if (!visited.insert(f).second)
          continue;
        for (int eId : sides[f])
        {
          auto& edge = e[eId];
          if (edge.a != static_cast<int>(v) && edge.b != static_cast<int>(v))
            continue;
          for (int neighbor : edge.faces)
            todo.push(neighbor);
        }
      }
      if (visited.size() != faces[v].size())
        fail(op, invalid_input, "nonmanifold disconnected vertex fan");
    }
  }
};

void validateCage(const Mesh& m)
{
  BRepBuilderAPI_Sewing sew(BRepExt::tolerance::length);
  for (const auto& q : m.q)
  {
    for (int side = 0; side < 4; ++side)
      if ((m.p[q[side]] - m.p[q[(side + 1) % 4]]).Modulus() < BRepExt::tolerance::length)
        fail(op, invalid_input, "degenerate cage edge (coincident or too close vertices)");
    for (const auto& tri :
         std::array<std::array<int, 3>, 2>{{{q[0], q[1], q[2]}, {q[0], q[2], q[3]}}})
    {
      if ((m.p[tri[1]] - m.p[tri[0]]).Crossed(m.p[tri[2]] - m.p[tri[0]]).Modulus() < 1e-14)
        fail(op, invalid_input, "degenerate cage quad triangle");
      BRepBuilderAPI_MakePolygon polygon;
      for (int v : tri)
        polygon.Add(gp_Pnt(m.p[v]));
      polygon.Close();
      BRepBuilderAPI_MakeFace face(polygon.Wire());
      if (!face.IsDone())
        fail(op, invalid_input, "invalid cage polygon");
      sew.Add(face.Face());
    }
  }
  sew.Perform();
  BRepAlgoAPI_Check check(sew.SewedShape(), false, true);
  if (!check.IsValid())
    fail(op, invalid_input, "the cage polygon surface self-intersects or is invalid");
}

Mesh refine(const Mesh& m)
{
  Topology t(m);
  Mesh     n;
  n.p.resize(m.p.size() + t.e.size() + m.q.size());
  n.corner.resize(n.p.size(), false);
  std::vector<gp_XYZ> centers;
  for (const auto& q : m.q)
    centers.push_back((m.p[q[0]] + m.p[q[1]] + m.p[q[2]] + m.p[q[3]]) / 4);
  for (size_t v = 0; v < m.p.size(); ++v)
  {
    n.corner[v] = m.corner[v];
    if (m.corner[v])
      n.p[v] = m.p[v];
    else if (!t.boundary[v].empty())
      n.p[v] = m.p[v] * 0.75 + (m.p[t.boundary[v][0]] + m.p[t.boundary[v][1]]) * 0.125;
    else
    {
      const double valence = t.faces[v].size();
      gp_XYZ       f(0, 0, 0), r(0, 0, 0);
      for (int id : t.faces[v])
        f += centers[id];
      for (int id : t.edges[v])
      {
        auto& e = t.e[id];
        r += (m.p[e.a] + m.p[e.b]) / 2;
      }
      n.p[v] = (f / valence + r * (2 / valence) + m.p[v] * (valence - 3)) / valence;
    }
  }
  for (size_t i = 0; i < t.e.size(); ++i)
  {
    const auto& e = t.e[i];
    gp_XYZ      p = m.p[e.a] + m.p[e.b];
    if (e.faces.size() == 2)
      p = (p + centers[e.faces[0]] + centers[e.faces[1]]) / 4;
    else
      p /= 2;
    n.p[m.p.size() + i] = p;
  }
  for (size_t f = 0; f < m.q.size(); ++f)
  {
    const auto& q = m.q[f];
    auto        e = t.sides[f];
    for (int& i : e)
      i += m.p.size();
    int c  = static_cast<int>(m.p.size() + t.e.size() + f);
    n.p[c] = centers[f];
    // All children keep the parent's UV orientation.
    n.q.push_back({q[0], e[0], c, e[3]});
    n.q.push_back({e[0], q[1], e[1], c});
    n.q.push_back({c, e[1], q[2], e[2]});
    n.q.push_back({e[3], c, e[2], q[3]});
  }
  return n;
}

std::vector<gp_XYZ> limitVertices(const Mesh& m)
{
  Topology            t(m);
  std::vector<gp_XYZ> p(m.p.size());
  for (size_t v = 0; v < m.p.size(); ++v)
  {
    if (m.corner[v])
      p[v] = m.p[v];
    else if (!t.boundary[v].empty())
      p[v] = (m.p[v] * 4 + m.p[t.boundary[v][0]] + m.p[t.boundary[v][1]]) / 6;
    else
    {
      const double n   = t.faces[v].size();
      gp_XYZ       sum = m.p[v] * (n * n);
      for (int eId : t.edges[v])
      {
        auto& e = t.e[eId];
        sum += m.p[e.a == static_cast<int>(v) ? e.b : e.a] * 4;
      }
      for (int f : t.faces[v])
      {
        const auto& q = m.q[f];
        int         i = 0;
        while (q[i] != static_cast<int>(v))
          ++i;
        sum += m.p[q[(i + 2) % 4]];
      }
      p[v] = sum / (n * (n + 5));
    }
  }
  return p;
}

// The dyadic UV lattice of each parent face; child orientation above makes it
// unambiguous.
std::vector<std::vector<gp_XYZ>> lattices(const Mesh& input, int levels)
{
  Mesh m = input;
  for (int i = 0; i < levels; ++i)
    m = refine(m);
  auto                             p         = limitVertices(m);
  int                              divisions = 1 << levels, stride = divisions + 1;
  std::vector<std::vector<gp_XYZ>> out(input.q.size(), std::vector<gp_XYZ>(stride * stride));
  for (size_t f = 0; f < m.q.size(); ++f)
  {
    size_t parent = f;
    int    u = 0, v = 0, step = 1;
    for (int l = 0; l < levels; ++l)
    {
      int child = parent % 4;
      parent /= 4;
      u += (child == 1 || child == 2) ? step : 0;
      v += child >= 2 ? step : 0;
      step *= 2;
    }
    auto&       grid               = out[parent];
    const auto& q                  = m.q[f];
    grid[u * stride + v]           = p[q[0]];
    grid[(u + 1) * stride + v]     = p[q[1]];
    grid[(u + 1) * stride + v + 1] = p[q[2]];
    grid[u * stride + v + 1]       = p[q[3]];
  }
  return out;
}

// Develop the oriented quad graph in integer UV. Holes, wrapping and
// extraordinary vertices cannot masquerade as one rectangular grid; each
// assignment is checked against all faces.
struct Grid
{
  int              nu = 0, nv = 0;
  std::vector<int> vertex, face;
};

Grid regularGrid(const Mesh& m, const Topology& t)
{
  std::map<int, Key> at;
  std::map<Key, int> occupied;
  std::queue<int>    todo;
  std::set<int>      seen;
  const int          initialX[4] = {0, 1, 1, 0}, initialY[4] = {0, 0, 1, 1};
  for (int i = 0; i < 4; ++i)
    at[m.q[0][i]] = {initialX[i], initialY[i]};
  todo.push(0);
  while (!todo.empty())
  {
    int f = todo.front();
    todo.pop();
    if (!seen.insert(f).second)
      continue;
    const auto& q = m.q[f];
    for (int i = 0; i < 4; ++i)
    {
      auto a = at[q[i]], b = at[q[(i + 1) % 4]], c = at[q[(i + 2) % 4]], d = at[q[(i + 3) % 4]];
      int  dx = b.first - a.first, dy = b.second - a.second;
      if (std::abs(dx) + std::abs(dy) != 1 || c != Key{b.first - dy, b.second + dx}
          || d != Key{a.first - dy, a.second + dx})
        return {};
      const auto& edge = t.e[t.sides[f][i]];
      if (edge.faces.size() != 2)
        continue;
      int         next = edge.faces[0] == f ? edge.faces[1] : edge.faces[0];
      const auto& n    = m.q[next];
      int         j    = 0;
      while (n[j] != q[(i + 1) % 4])
        ++j;
      std::array<Key, 4> coords = {b,
                                   a,
                                   Key{a.first + dy, a.second - dx},
                                   Key{b.first + dy, b.second - dx}};
      for (int k = 0; k < 4; ++k)
      {
        int  id    = n[(j + k) % 4];
        auto found = at.emplace(id, coords[k]);
        if (!found.second && found.first->second != coords[k])
          return {};
      }
      todo.push(next);
    }
  }
  if (at.size() != m.p.size() || seen.size() != m.q.size())
    return {};
  int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
  for (auto [id, xy] : at)
  {
    if (!occupied.emplace(xy, id).second)
      return {};
    x0 = std::min(x0, xy.first);
    x1 = std::max(x1, xy.first);
    y0 = std::min(y0, xy.second);
    y1 = std::max(y1, xy.second);
  }
  Grid g;
  g.nu = x1 - x0;
  g.nv = y1 - y0;
  if (g.nu * g.nv != static_cast<int>(m.q.size())
      || (g.nu + 1) * (g.nv + 1) != static_cast<int>(m.p.size()))
    return {};
  g.vertex.resize(m.p.size());
  g.face.resize(m.q.size());
  for (auto [xy, id] : occupied)
    g.vertex[(xy.first - x0) * (g.nv + 1) + xy.second - y0] = id;
  for (size_t f = 0; f < m.q.size(); ++f)
  {
    int x = x1, y = y1;
    for (int id : m.q[f])
    {
      x = std::min(x, at[id].first);
      y = std::min(y, at[id].second);
    }
    g.face[(x - x0) * g.nv + y - y0] = static_cast<int>(f);
  }
  return g;
}

SplineNet bezierNet()
{
  SplineNet n;
  n.degreeU = n.degreeV = 3;
  n.knotsU = n.knotsV = {0, 0, 0, 0, 1, 1, 1, 1};
  n.poles.resize(16);
  return n;
}

gp_XYZ gridPole(const Mesh& m, const Grid& g, int u, int v)
{
  if (u < 0)
    return gridPole(m, g, 0, v) * 2 - gridPole(m, g, 1, v);
  if (u > g.nu)
    return gridPole(m, g, g.nu, v) * 2 - gridPole(m, g, g.nu - 1, v);
  if (v < 0)
    return gridPole(m, g, u, 0) * 2 - gridPole(m, g, u, 1);
  if (v > g.nv)
    return gridPole(m, g, u, g.nv) * 2 - gridPole(m, g, u, g.nv - 1);
  return m.p[g.vertex[u * (g.nv + 1) + v]];
}

std::vector<SplineNet> exactGrid(const Mesh& m, const Grid& g)
{
  const double           b[4][4] = {{1. / 6, 4. / 6, 1. / 6, 0},
                                    {0, 4. / 6, 2. / 6, 0},
                                    {0, 2. / 6, 4. / 6, 0},
                                    {0, 1. / 6, 4. / 6, 1. / 6}};
  std::vector<SplineNet> out(m.q.size());
  for (int u = 0; u < g.nu; ++u)
    for (int v = 0; v < g.nv; ++v)
    {
      auto n = bezierNet();
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
          gp_XYZ p(0, 0, 0);
          for (int k = 0; k < 4; ++k)
            for (int l = 0; l < 4; ++l)
              p += gridPole(m, g, u + k - 1, v + l - 1) * (b[i][k] * b[j][l]);
          n.poles[i * 4 + j] = p;
        }
      // Developed UV may be rotated from this face's input UV.
      int         f      = g.face[u * g.nv + v];
      const auto& q      = m.q[f];
      const int   ids[4] = {g.vertex[u * (g.nv + 1) + v],
                            g.vertex[(u + 1) * (g.nv + 1) + v],
                            g.vertex[(u + 1) * (g.nv + 1) + v + 1],
                            g.vertex[u * (g.nv + 1) + v + 1]};
      int         r      = 0;
      while (q[0] != ids[r])
        ++r;
      auto rotated = n;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
          int x = i, y = j;
          if (r == 1)
          {
            x = 3 - j;
            y = i;
          }
          if (r == 2)
          {
            x = 3 - i;
            y = 3 - j;
          }
          if (r == 3)
          {
            x = j;
            y = 3 - i;
          }
          rotated.poles[i * 4 + j] = n.poles[x * 4 + y];
        }
      out[f] = rotated;
    }
  return out;
}

using Curve = std::array<gp_XYZ, 4>;

std::vector<SplineNet> fitted(const Mesh&                             m,
                              const Topology&                         t,
                              const std::vector<std::vector<gp_XYZ>>& samples)
{
  std::vector<Curve> curves(t.e.size());
  for (size_t e = 0; e < t.e.size(); ++e)
  {
    const auto& edge = t.e[e];
    int         f = edge.faces[0], side = 0;
    while (t.sides[f][side] != static_cast<int>(e))
      ++side;
    auto along = [&](int k) {
      int u = 0, v = 0;
      if (side == 0)
        u = k;
      if (side == 1)
      {
        u = 4;
        v = k;
      }
      if (side == 2)
      {
        u = 4 - k;
        v = 4;
      }
      if (side == 3)
        v = 4 - k;
      return samples[f][u * 5 + v];
    };
    Curve c;
    c[0]      = along(0);
    c[3]      = along(4);
    double aa = 0, ab = 0;
    gp_XYZ a(0, 0, 0), b(0, 0, 0);
    for (int k = 1; k < 4; ++k)
    {
      double x = k / 4., y = 1 - x, b1 = 3 * x * y * y, b2 = 3 * x * x * y;
      auto   rhs = along(k) - c[0] * (y * y * y) - c[3] * (x * x * x);
      aa += b1 * b1;
      ab += b1 * b2;
      a += rhs * b1;
      b += rhs * b2;
    }
    double determinant = aa * aa - ab * ab;
    c[1]               = (a * aa - b * ab) / determinant;
    c[2]               = (b * aa - a * ab) / determinant;
    if (m.q[f][side] != edge.a)
      std::reverse(c.begin(), c.end());
    curves[e] = c;
  }
  std::vector<SplineNet> out;
  for (size_t f = 0; f < m.q.size(); ++f)
  {
    auto              n = bezierNet();
    std::vector<char> fixed(16, false);
    for (int side = 0; side < 4; ++side)
    {
      auto c = curves[t.sides[f][side]];
      if (m.q[f][side] != t.e[t.sides[f][side]].a)
        std::reverse(c.begin(), c.end());
      for (int k = 0; k < 4; ++k)
      {
        int u = 0, v = 0;
        if (side == 0)
          u = k;
        if (side == 1)
        {
          u = 3;
          v = k;
        }
        if (side == 2)
        {
          u = 3 - k;
          v = 3;
        }
        if (side == 3)
          v = 3 - k;
        n.poles[u * 4 + v] = c[k];
        fixed[u * 4 + v]   = true;
      }
    }
    std::vector<FitCondition> conditions;
    for (int i = 1; i < 4; ++i)
      for (int j = 1; j < 4; ++j)
        for (int axis = 0; axis < 3; ++axis)
        {
          gp_XYZ d(0, 0, 0);
          d.SetCoord(axis + 1, 1);
          FitCondition c{};
          c.u        = i / 4.;
          c.v        = j / 4.;
          c.count    = 1;
          c.terms[0] = {0, 0, d};
          c.value    = samples[f][i * 5 + j].Coord(axis + 1);
          c.weight   = 1;
          c.soft     = true;
          conditions.push_back(c);
        }
    // The single shared fitter solves only the four free interior poles. Every
    // boundary coefficient is fixed to the common edge curve, so independent
    // fits cannot crack.
    fitNet(op, n, fixed, conditions, {0, 0}, {});
    out.push_back(n);
  }
  return out;
}

gp_XYZ normal(const occ::handle<Geom_BSplineSurface>& s, double u, double v, double* jac = nullptr)
{
  gp_Pnt p;
  gp_Vec du, dv;
  s->D1(u, v, p, du, dv);
  gp_XYZ n      = du.XYZ().Crossed(dv.XYZ());
  double length = n.Modulus();
  if (jac)
    *jac = length;
  if (!std::isfinite(length) || length < 1e-14)
    fail(op, no_solution, "singular patch derivative");
  return n / length;
}

std::pair<double, double> edgeUV(int side, double t)
{
  if (side == 0)
    return {t, 0};
  if (side == 1)
    return {1, t};
  if (side == 2)
    return {1 - t, 1};
  return {0, 1 - t};
}

occ::handle<Geom_BSplineSurface> mergeGrid(const std::vector<SplineNet>& nets,
                                           const Mesh&                   m,
                                           const Grid&                   g)
{
  // Assemble exactly matching Bezier coefficients, then remove only the regular
  // knots.
  int                        nu = 3 * g.nu + 1, nv = 3 * g.nv + 1;
  NCollection_Array2<gp_Pnt> p(1, nu, 1, nv);
  for (int u = 0; u < g.nu; ++u)
    for (int v = 0; v < g.nv; ++v)
    {
      int         f     = g.face[u * g.nv + v];
      const auto& q     = m.q[f];
      int         first = g.vertex[u * (g.nv + 1) + v];
      int         r     = 0;
      while (q[r] != first)
        ++r;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
          int x = i, y = j;
          if (r == 1)
          {
            x = 3 - j;
            y = i;
          }
          if (r == 2)
          {
            x = 3 - i;
            y = 3 - j;
          }
          if (r == 3)
          {
            x = j;
            y = 3 - i;
          }
          p.SetValue(3 * u + i + 1, 3 * v + j + 1, gp_Pnt(nets[f].poles[x * 4 + y]));
        }
    }
  NCollection_Array1<double> ku(1, g.nu + 1), kv(1, g.nv + 1);
  NCollection_Array1<int>    mu(1, g.nu + 1), mv(1, g.nv + 1);
  for (int i = 1; i <= g.nu + 1; ++i)
  {
    ku.SetValue(i, (i - 1.) / g.nu);
    mu.SetValue(i, i == 1 || i == g.nu + 1 ? 4 : 3);
  }
  for (int i = 1; i <= g.nv + 1; ++i)
  {
    kv.SetValue(i, (i - 1.) / g.nv);
    mv.SetValue(i, i == 1 || i == g.nv + 1 ? 4 : 3);
  }
  auto   s     = occ::handle<Geom_BSplineSurface>(new Geom_BSplineSurface(p, ku, kv, mu, mv, 3, 3));
  double scale = 0;
  for (const auto& xyz : m.p)
    scale = std::max(scale, xyz.Modulus());
  double epsilon = std::max(1e-12, scale * 1e-12);
  for (int i = 2; i < s->NbUKnots(); ++i)
    if (!s->RemoveUKnot(i, 1, epsilon))
      fail(op, no_solution, "regular U patch merge failed");
  for (int i = 2; i < s->NbVKnots(); ++i)
    if (!s->RemoveVKnot(i, 1, epsilon))
      fail(op, no_solution, "regular V patch merge failed");
  return s;
}
} // namespace

PolySplineResult polySpline(OcctKernel&                kernel,
                            const std::vector<double>& xyz,
                            const std::vector<int>&    quads,
                            double                     tolerance,
                            int                        refinement,
                            bool                       merge,
                            bool                       rounded,
                            bool                       boundaryExactly,
                            bool                       requireG1,
                            double                     angularTolerance)
{
  if (xyz.size() < 12 || xyz.size() % 3 || xyz.size() > 1024 * 3 || quads.empty()
      || quads.size() % 4 || quads.size() > 256 * 4)
    fail(op, invalid_input, "cages need XYZ vertices and quads, at most 1024 vertices / 256 faces");
  if (!std::isfinite(tolerance) || tolerance < BRepExt::tolerance::length || tolerance > 1
      || refinement < 0 || refinement > 3 || !std::isfinite(angularTolerance)
      || angularTolerance < 0 || angularTolerance > 90)
    fail(op, invalid_input, "invalid tolerance, angular tolerance or refinement budget");
  if (rounded)
    fail(op,
         unsupported_geometry,
         "Rounded Corners is not implemented; the V1 boundary keeps one-face "
         "corners fixed");
  Mesh cage;
  for (size_t i = 0; i < xyz.size(); i += 3)
  {
    for (int a = 0; a < 3; ++a)
      if (!std::isfinite(xyz[i + a]) || std::abs(xyz[i + a]) > 1e9)
        fail(op, invalid_input, "cage coordinates must be finite and bounded");
    cage.p.emplace_back(xyz[i], xyz[i + 1], xyz[i + 2]);
  }
  for (size_t i = 0; i < quads.size(); i += 4)
  {
    Quad q{};
    for (int a = 0; a < 4; ++a)
    {
      q[a] = quads[i + a];
      if (q[a] < 0 || q[a] >= static_cast<int>(cage.p.size()))
        fail(op, invalid_input, "quad names an absent vertex");
    }
    cage.q.push_back(q);
  }
  Topology original(cage);
  cage.corner.resize(cage.p.size(), false);
  PolySplineResult report;
  std::set<int>    extraordinary;
  for (size_t v = 0; v < cage.p.size(); ++v)
  {
    if (!original.boundary[v].empty())
    {
      if (original.faces[v].size() > 2)
        fail(op, unsupported_geometry, "extraordinary boundary vertices are not supported");
      cage.corner[v] = original.faces[v].size() == 1;
    }
    else if (original.faces[v].size() != 4)
    {
      int n = original.faces[v].size();
      if (n != 3 && n != 5 && n != 6)
        fail(op, unsupported_geometry, "only isolated interior valences 3, 5 and 6 are supported");
      extraordinary.insert(v);
    }
  }
  for (int v : extraordinary)
    for (int eId : original.edges[v])
    {
      const auto& e     = original.e[eId];
      int         other = e.a == v ? e.b : e.a;
      if (extraordinary.count(other))
        fail(op,
             unsupported_geometry,
             "adjacent extraordinary vertices are not supported; add a genuine "
             "geometric edge loop");
      if (!original.boundary[other].empty())
        fail(op,
             unsupported_geometry,
             "extraordinary vertices next to a boundary are not supported");
    }
  std::set<int>   visited;
  std::queue<int> todo;
  todo.push(0);
  while (!todo.empty())
  {
    int f = todo.front();
    todo.pop();
    if (!visited.insert(f).second)
      continue;
    for (int e : original.sides[f])
      for (int n : original.e[e].faces)
        todo.push(n);
  }
  if (visited.size() != cage.q.size())
    fail(op,
         unsupported_geometry,
         "disconnected components must be imported as separate OBJ objects in "
         "V1");
  validateCage(cage);
  report.extraordinary = extraordinary.size();
  auto grid            = regularGrid(cage, original);
  report.regular       = grid.nu > 0;
  if (merge && !report.regular)
    fail(
      op,
      unsupported_geometry,
      "Merge Patches is supported only for rectangular regular grids; disable it for this cage.");
  if (boundaryExactly && !report.regular)
    fail(op,
         unsupported_geometry,
         "Boundary Exactly is supported only for a regular rectangular grid in "
         "V1");
  Mesh                   work = cage;
  std::vector<SplineNet> nets;
  Grid                   workingGrid = grid;
  if (report.regular)
  {
    nets       = exactGrid(cage, grid);
    refinement = 0;
  }
  else
  {
    int factor = 1 << (2 * refinement);
    if (cage.q.size() * factor > 1024)
      fail(op,
           unsupported_geometry,
           "conversion exceeds the 1024-patch budget; reduce refinement or "
           "cage size");
    for (int i = 0; i < refinement; ++i)
      work = refine(work);
    Topology t(work);
    nets = fitted(work, t, lattices(work, 2));
  }
  report.refinement = refinement;
  Topology                                      topology(work);
  auto                                          oracle = lattices(work, 3);
  std::vector<occ::handle<Geom_BSplineSurface>> surfaces;
  for (const auto& net : nets)
    surfaces.push_back(surfaceOf(net));
  report.minJacobian = std::numeric_limits<double>::infinity();
  for (size_t f = 0; f < surfaces.size(); ++f)
  {
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; j <= 8; ++j)
      {
        auto   s         = surfaces[f];
        auto   p         = s->Value(i / 8., j / 8.).XYZ();
        double gap       = (p - oracle[f][i * 9 + j]).Modulus();
        report.deviation = std::max(report.deviation, gap);
        ++report.samples;
        // The extraordinary corner is parametrically singular; test open
        // interior and seams away from endpoints. Its normal is not invented
        // from shaded triangles.
        if (i > 0 && i < 8 && j > 0 && j < 8)
        {
          double jac;
          auto   n           = normal(s, i / 8., j / 8., &jac);
          report.minJacobian = std::min(report.minJacobian, jac);
          auto u             = oracle[f][(i + 1) * 9 + j] - oracle[f][(i - 1) * 9 + j],
               v             = oracle[f][i * 9 + j + 1] - oracle[f][i * 9 + j - 1];
          if (n.Dot(u.Crossed(v)) <= 0)
            fail(op,
                 no_solution,
                 "patch folds or reverses relative to the subdivision limit "
                 "samples");
        }
      }
  }
  if (report.deviation > tolerance)
  {
    std::ostringstream msg;
    msg << "sampled limit deviation " << report.deviation << " mm exceeds " << tolerance
        << " mm; increase bounded refinement";
    fail(op, no_solution, msg.str());
  }
  for (size_t e = 0; e < topology.e.size(); ++e)
  {
    const auto& edge = topology.e[e];
    if (edge.faces.size() != 2)
      continue;
    ++report.seams;
    int a = edge.faces[0], b = edge.faces[1], sa = 0, sb = 0;
    while (topology.sides[a][sa] != static_cast<int>(e))
      ++sa;
    while (topology.sides[b][sb] != static_cast<int>(e))
      ++sb;
    for (int k = 0; k <= 32; ++k)
    {
      double t  = k / 32.;
      auto   uv = edgeUV(sa, t), xy = edgeUV(sb, 1 - t);
      auto p = surfaces[a]->Value(uv.first, uv.second), q = surfaces[b]->Value(xy.first, xy.second);
      report.gap = std::max(report.gap, p.Distance(q));
      if (k > 0 && k < 32)
      {
        auto n             = normal(surfaces[a], uv.first, uv.second),
             m             = normal(surfaces[b], xy.first, xy.second);
        double angle       = std::atan2(n.Crossed(m).Modulus(), n.Dot(m)) * 180 / std::acos(-1.);
        report.normalAngle = std::max(report.normalAngle, angle);
      }
    }
  }
  if (report.gap > std::min(tolerance, 1e-6))
    fail(op, no_solution, "shared patch edge gap exceeds the seam tolerance");
  if (requireG1 && report.normalAngle > angularTolerance)
  {
    std::ostringstream msg;
    msg << "measured seam normal angle " << report.normalAngle << " degrees exceeds requested "
        << angularTolerance << " degrees; G1 was refused";
    fail(op, no_solution, msg.str());
  }
  BRepBuilderAPI_Sewing    sew(BRepExt::tolerance::length);
  std::vector<TopoDS_Face> faces;
  if (merge && report.regular)
  {
    auto s = mergeGrid(nets, cage, grid);
    for (int u = 0; u < grid.nu; ++u)
      for (int v = 0; v < grid.nv; ++v)
        for (int i = 0; i <= 8; ++i)
          for (int j = 0; j <= 8; ++j)
          {
            auto p = s->Value((u + i / 8.) / grid.nu, (v + j / 8.) / grid.nv);
            // Compare the merged representation against the exact pre-merge
            // patches.
            int f = grid.face[u * grid.nv + v];
            int r = 0;
            while (cage.q[f][r] != grid.vertex[u * (grid.nv + 1) + v])
              ++r;
            double x = i / 8., y = j / 8.;
            if (r == 1)
            {
              double old = x;
              x          = 1 - y;
              y          = old;
            }
            if (r == 2)
            {
              x = 1 - x;
              y = 1 - y;
            }
            if (r == 3)
            {
              double old = x;
              x          = y;
              y          = 1 - old;
            }
            if (p.Distance(surfaces[f]->Value(x, y)) > 1e-8)
              fail(op, no_solution, "merged patch representation changed the limit surface");
          }
    BRepBuilderAPI_MakeFace make(s, BRepExt::tolerance::length);
    if (!make.IsDone())
      fail(op, no_solution, "merged face construction failed");
    faces.push_back(make.Face());
    report.provenance.push_back(-1);
    // A merged multispan face has no B-rep joins. The pre-merge regular
    // joins were checked above; do not report them as output seams.
    report.seams       = 0;
    report.gap         = 0;
    report.normalAngle = 0;
  }
  else
  {
    for (size_t f = 0; f < surfaces.size(); ++f)
    {
      BRepBuilderAPI_MakeFace make(surfaces[f], BRepExt::tolerance::length);
      if (!make.IsDone())
        fail(op, no_solution, "patch face construction failed");
      faces.push_back(make.Face());
      report.provenance.push_back(static_cast<int>(f) / (1 << (2 * refinement)));
    }
  }
  for (const auto& f : faces)
    sew.Add(f);
  sew.Perform();
  if (sew.NbMultipleEdges())
    fail(op, no_solution, "sewing found multiple edges");
  report.freeEdges    = sew.NbFreeEdges();
  TopoDS_Shape result = sew.SewedShape();
  bool         closed = std::all_of(original.e.begin(), original.e.end(), [](const Edge& e) {
    return e.faces.size() == 2;
  });
  if (closed)
  {
    if (report.freeEdges || result.ShapeType() != TopAbs_SHELL)
      fail(op, no_solution, "closed cage did not produce one closed shell");
    BRepBuilderAPI_MakeSolid solid(TopoDS::Shell(result));
    result        = solid.Solid();
    auto oriented = TopoDS::Solid(result);
    if (!BRepLib::OrientClosedSolid(oriented))
      fail(op, no_solution, "closed shell cannot be oriented as a solid");
    result = oriented;
    GProp_GProps props;
    BRepGProp::VolumeProperties(result, props);
    if (props.Mass() <= 1e-12)
      fail(op, no_solution, "closed shell has no positive volume");
    report.solid = true;
  }
  if (result.IsNull())
    fail(op, no_solution, "generated B-rep is empty");
  BRepCheck_Analyzer analyzer(result, true);
  analyzer.SetExactMethod(true);
  if (!analyzer.IsValid())
    fail(op, no_solution, "generated B-rep is invalid (including edge/pcurve agreement)");
  const auto sources = report.provenance;
  report.provenance.clear();
  for (TopExp_Explorer it(result, TopAbs_FACE); it.More(); it.Next())
  {
    bool found = false;
    for (size_t f = 0; f < faces.size(); ++f)
    {
      auto product = sew.IsModified(faces[f]) ? sew.Modified(faces[f]) : faces[f];
      if (it.Current().IsSame(product))
      {
        report.provenance.push_back(sources[f]);
        found = true;
        break;
      }
    }
    if (!found)
      fail(op, no_solution, "sewing lost patch-to-cage provenance");
  }
  BRepAlgoAPI_Check check(result, false, true);
  if (!check.IsValid())
    fail(op, no_solution, "generated B-rep has self-interference or invalid geometry");
  report.patches  = faces.size();
  report.resultId = store(kernel, result);
  return report;
}
} // namespace BRepExt
