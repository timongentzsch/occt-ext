// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepAlgoExt_PatchFlow.hxx>
#include <GeomAlgoExt_Phase.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <GeomAlgoExt_Resolution.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
#include <BRepMesh_Delaun.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_Surface.hxx>
#include <NCollection_IncAllocator.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

namespace loft {
namespace {
constexpr const char *op = "patch";

/// The most nodes the triangulation is given. The flattening's linear solves cost about this
/// to the power 1.5; a chart finer than this is not what limits the fit.
constexpr int mostNodes = 32000;

/// How near the conjugate gradients bring a solve, as a share of its right-hand side.
constexpr double solved = 1e-13;

/// The cotangent Laplacian of a triangulation in space: per node its neighbours and weights,
/// and the sum of its weights.
struct Laplacian {
    std::vector<std::vector<std::pair<int, double>>> neighbours;
    std::vector<double> diagonal;
    /// (A x) at one node.
    double row(const std::vector<double> &x, int i) const {
        double sum = diagonal[i] * x[i];
        for (const auto &[j, w] : neighbours[i]) sum -= w * x[j];
        return sum;
    }
};

/// Solve A x = b for the nodes that are not `fixed`, with x as given where they are, by
/// conjugate gradients preconditioned by the diagonal.
void solve(const Laplacian &a, const std::vector<char> &fixed, const std::vector<double> &b, std::vector<double> &x) {
    const size_t n = x.size();
    std::vector<double> r(n, 0), z(n, 0), p(n, 0), q(n, 0);
    double norm = 0;
    for (size_t i = 0; i < n; ++i)
        if (!fixed[i]) {
            r[i] = b[i] - a.row(x, static_cast<int>(i));
            norm += b[i] * b[i];
        }
    // Where the right-hand side is nothing, the fixed values drive the solve.
    double start = 0;
    for (size_t i = 0; i < n; ++i) start += r[i] * r[i];
    const double target = solved * solved * std::max(norm, start);
    if (!(start > target)) return;
    double rz = 0;
    for (size_t i = 0; i < n; ++i)
        if (!fixed[i]) {
            z[i] = r[i] / a.diagonal[i];
            p[i] = z[i];
            rz += r[i] * z[i];
        }
    // The free part of the matrix in one array, row after row: the passes read only it.
    std::vector<size_t> begin(n + 1, 0);
    std::vector<std::pair<int, double>> free;
    for (size_t i = 0; i < n; ++i) {
        begin[i] = free.size();
        if (!fixed[i])
            for (const auto &[j, w] : a.neighbours[i])
                if (!fixed[j]) free.push_back({j, w});
    }
    begin[n] = free.size();
    for (size_t pass = 0; pass < 4 * n; ++pass) {
        double pq = 0;
        for (size_t i = 0; i < n; ++i)
            if (!fixed[i]) {
                q[i] = a.diagonal[i] * p[i];
                for (size_t k = begin[i]; k < begin[i + 1]; ++k) q[i] -= free[k].second * p[free[k].first];
                pq += p[i] * q[i];
            }
        if (!(pq > 0)) fail(op, no_solution, "the conformal chart's equations are not positive definite");
        const double step = rz / pq;
        double left = 0, next = 0;
        for (size_t i = 0; i < n; ++i)
            if (!fixed[i]) {
                x[i] += step * p[i];
                r[i] -= step * q[i];
                left += r[i] * r[i];
                z[i] = r[i] / a.diagonal[i];
                next += r[i] * z[i];
            }
        if (left <= target) return;
        for (size_t i = 0; i < n; ++i)
            if (!fixed[i]) p[i] = z[i] + next / rz * p[i];
        rz = next;
    }
    fail(op, no_solution, "the conformal chart's equations do not converge");
}

double signedArea(const gp_XY &a, const gp_XY &b, const gp_XY &c) { return (b - a).Crossed(c - a) / 2; }

} // namespace

ConformalChart conformalChart(const TopoDS_Face &face, double spacing, double deflection) {
    PatchPhase phase("conformal-chart");
    ConformalChart chart;
    const occ::handle<Geom_Surface> surface = BRep_Tool::Surface(face);
    // The boundary in the face's parameters: each edge at even steps of its length in space.
    std::vector<gp_XY> rim;
    const auto divide = [&](double step) {
        rim.clear();
        for (BRepTools_WireExplorer it(BRepTools::OuterWire(face), face); it.More(); it.Next()) {
            const TopoDS_Edge &edge = it.Current();
            double first, last;
            const auto curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
            if (curve.IsNull()) fail(op, no_solution, "the face to be flattened has an edge without a curve on it");
            BRepAdaptor_Curve inSpace(edge);
            const double length = GCPnts_AbscissaPoint::Length(inSpace, first, last);
            const int steps = std::max(2, static_cast<int>(std::ceil(length / step)));
            GCPnts_UniformAbscissa even(inSpace, steps + 1, first, last);
            if (!even.IsDone() || even.NbPoints() != steps + 1)
                fail(op, no_solution, "the face to be flattened has an edge that cannot be divided");
            const bool reversed = edge.Orientation() == TopAbs_REVERSED;
            for (int k = 0; k < steps; ++k)
                rim.push_back(curve->Value(even.Parameter(reversed ? steps + 1 - k : 1 + k)).XY());
        }
        double turning = 0;
        for (size_t i = 0; i < rim.size(); ++i) turning += rim[i].Crossed(rim[(i + 1) % rim.size()]);
        if (turning < 0) std::reverse(rim.begin(), rim.end());
        return std::abs(turning) / 2;
    };
    // Whether face parameters are inside the boundary, and further from it than `margin`.
    const auto inside = [&](const gp_XY &p, double margin) {
        bool in = false;
        double nearest = HUGE_VAL;
        for (size_t i = 0, j = rim.size() - 1; i < rim.size(); j = i++) {
            const gp_XY edge = rim[i] - rim[j], to = p - rim[j];
            nearest = std::min(nearest, (to - edge * std::clamp(to.Dot(edge) / edge.SquareModulus(), 0.0, 1.0)).Modulus());
            if ((rim[i].Y() > p.Y()) != (rim[j].Y() > p.Y()) && p.X() < rim[j].X() + edge.X() * to.Y() / edge.Y())
                in = !in;
        }
        return in && nearest > margin;
    };
    const auto bounds = [&] {
        std::pair<gp_XY, gp_XY> out{rim.front(), rim.front()};
        for (const gp_XY &p : rim) {
            out.first.SetCoord(std::min(out.first.X(), p.X()), std::min(out.first.Y(), p.Y()));
            out.second.SetCoord(std::max(out.second.X(), p.X()), std::max(out.second.Y(), p.Y()));
        }
        return out;
    };
    // Where the face bends inside the boundary, triangles short enough to stay within
    // `deflection` of it (a chord of length l on a bend k stands k l^2 / 8 off), as far as the
    // most nodes allow.
    const double areaInParameters = divide(spacing);
    double bend = 0;
    {
        const auto [low, high] = bounds();
        for (double u = low.X(); u <= high.X(); u += spacing)
            for (double v = low.Y(); v <= high.Y(); v += spacing) {
                if (!inside(gp_XY(u, v), 0)) continue;
                GeomLProp_SLProps props(surface, u, v, 2, loft::tolerance::parameterMerge);
                if (props.IsCurvatureDefined())
                    bend = std::max({bend, std::abs(props.MaxCurvature()), std::abs(props.MinCurvature())});
            }
    }
    if (bend > 0) spacing = std::min(spacing, std::sqrt(8 * deflection / bend));
    spacing = std::max(spacing, std::sqrt(areaInParameters / (mostNodes * std::sqrt(3.0) / 4)));
    divide(spacing);
    // Inside, an even triangular lattice, kept half a step off the boundary.
    const auto [from, to] = bounds();
    std::vector<gp_XY> nodes = rim;
    const double rowStep = spacing * std::sqrt(3.0) / 2;
    for (int row = 0; from.Y() + row * rowStep < to.Y(); ++row)
        for (double u = from.X() + (row % 2) * spacing / 2; u < to.X(); u += spacing)
            if (const gp_XY p(u, from.Y() + row * rowStep); inside(p, spacing / 2)) nodes.push_back(p);
    if (static_cast<int>(nodes.size()) > mostNodes)
        fail(op, no_solution, "the face to be flattened needs more nodes than the chart is given");
    // Constrained Delaunay in the face's parameters, with the boundary as its frontier.
    occ::handle<NCollection_IncAllocator> allocator = new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
    occ::handle<BRepMesh_DataStructureOfDelaun> structure = new BRepMesh_DataStructureOfDelaun(allocator, static_cast<int>(nodes.size()));
    std::vector<int> ids;
    for (size_t i = 0; i < nodes.size(); ++i)
        ids.push_back(structure->AddNode(BRepMesh_Vertex(nodes[i], static_cast<int>(i), i < rim.size() ? BRepMesh_Frontier : BRepMesh_Free), true));
    for (size_t i = 0; i < rim.size(); ++i)
        structure->AddLink(BRepMesh_Edge(ids[i], ids[(i + 1) % rim.size()], BRepMesh_Frontier));
    IMeshData::VectorOfInteger order(static_cast<int>(nodes.size()), allocator);
    for (int id : ids) order.Append(id);
    BRepMesh_Delaun mesher(structure, order);
    std::unordered_map<int, int> node;
    for (size_t i = 0; i < ids.size(); ++i) node[ids[i]] = static_cast<int>(i);
    for (IMeshData::IteratorOfMapOfInteger it(structure->ElementsOfDomain()); it.More(); it.Next()) {
        int corner[3];
        structure->ElementNodes(structure->GetElement(it.Key()), corner);
        std::array<int, 3> t{};
        for (int k = 0; k < 3; ++k) {
            const auto found = node.find(corner[k]);
            if (found == node.end()) fail(op, no_solution, "the face's triangulation has a node of its own");
            t[k] = found->second;
        }
        chart.triangles.push_back(t);
    }
    // A free node that the triangulation leaves out, as it does one too near the boundary to
    // insert, is no part of the chart: only the nodes of its triangles are kept, in their order.
    std::vector<int> kept(nodes.size(), -1);
    for (const auto &t : chart.triangles)
        for (int k : t) kept[k] = 0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (kept[i] < 0) continue;
        kept[i] = static_cast<int>(chart.source.size());
        chart.source.push_back(nodes[i]);
        chart.point.push_back(surface->Value(nodes[i].X(), nodes[i].Y()));
    }
    for (auto &t : chart.triangles)
        for (int &k : t) k = kept[k];
    // Every triangle anticlockwise in the face's parameters, which is anticlockwise about the
    // surface's own normal.
    for (auto &t : chart.triangles) {
        const double area = signedArea(chart.source[t[0]], chart.source[t[1]], chart.source[t[2]]);
        if (!(std::abs(area) > 0)) fail(op, no_solution, "the face's triangulation has a degenerate triangle");
        if (area < 0) std::swap(t[1], t[2]);
    }
    const int n = static_cast<int>(chart.source.size());

    // The boundary: the edges that one triangle has and no other, in turn.
    std::unordered_map<long long, int> directed;
    for (const auto &t : chart.triangles)
        for (int k = 0; k < 3; ++k) directed[static_cast<long long>(t[k]) * n + t[(k + 1) % 3]] = t[(k + 2) % 3];
    std::vector<int> next(n, -1);
    int edges = 0;
    for (const auto &[key, opposite] : directed) {
        const int a = static_cast<int>(key / n), b = static_cast<int>(key % n);
        if (directed.count(static_cast<long long>(b) * n + a)) continue;
        if (next[a] >= 0) fail(op, no_solution, "the face's triangulation is not one disc");
        next[a] = b;
        ++edges;
    }
    if (!edges) fail(op, no_solution, "the face's triangulation has no boundary");
    const int first = static_cast<int>(std::find_if(next.begin(), next.end(), [](int k) { return k >= 0; }) - next.begin());
    for (int at = first;;) {
        chart.boundary.push_back(at);
        at = next[at];
        if (at < 0 || static_cast<int>(chart.boundary.size()) > edges)
            fail(op, no_solution, "the face's triangulation is not one disc");
        if (at == first) break;
    }
    if (static_cast<int>(chart.boundary.size()) != edges)
        fail(op, no_solution, "the face's triangulation is not one disc");
    std::vector<char> onBoundary(n, 0);
    for (int b : chart.boundary) onBoundary[b] = 1;

    // Cotangent weights and the angles at every node, in space.
    Laplacian laplacian{std::vector<std::vector<std::pair<int, double>>>(n), std::vector<double>(n, 0)};
    std::unordered_map<long long, double> weight;
    std::vector<double> angles(n, 0);
    for (const auto &t : chart.triangles)
        for (int k = 0; k < 3; ++k) {
            const int i = t[k], j = t[(k + 1) % 3], l = t[(k + 2) % 3];
            const gp_Vec a(chart.point[l], chart.point[i]), b(chart.point[l], chart.point[j]);
            const double cross = a.Crossed(b).Magnitude(), dot = a.Dot(b);
            if (!(cross > 0)) fail(op, no_solution, "the face's triangulation has a degenerate triangle in space");
            angles[l] += std::atan2(cross, dot);
            weight[static_cast<long long>(std::min(i, j)) * n + std::max(i, j)] += dot / cross / 2;
        }
    for (const auto &[key, w] : weight) {
        const int i = static_cast<int>(key / n), j = static_cast<int>(key % n);
        laplacian.neighbours[i].push_back({j, w});
        laplacian.neighbours[j].push_back({i, w});
        laplacian.diagonal[i] += w;
        laplacian.diagonal[j] += w;
    }
    for (int i = 0; i < n; ++i)
        if (!(laplacian.diagonal[i] > 0)) fail(op, no_solution, "the face's triangulation has a node of no weight");

    // The scale that takes the curvature out of the inside, nothing along the boundary; and
    // the boundary's turning that comes of it.
    std::vector<double> scale(n, 0), curvature(n, 0);
    for (int i = 0; i < n; ++i) curvature[i] = onBoundary[i] ? 0 : 2 * M_PI - angles[i];
    solve(laplacian, onBoundary, curvature, scale);
    const size_t m = chart.boundary.size();
    std::vector<double> length(m), heading(m);
    double turned = 0;
    for (size_t k = 0; k < m; ++k) {
        const int b = chart.boundary[k];
        length[k] = chart.point[b].Distance(chart.point[chart.boundary[(k + 1) % m]]);
        if (k > 0) turned += M_PI - angles[b] - laplacian.row(scale, b);
        heading[k] = turned;
    }
    // Close the boundary by the least change of its lengths, each counted by its own length.
    double gxx = 0, gxy = 0, gyy = 0, rx = 0, ry = 0;
    for (size_t k = 0; k < m; ++k) {
        const double c = std::cos(heading[k]), s = std::sin(heading[k]);
        gxx += length[k] * c * c; gxy += length[k] * c * s; gyy += length[k] * s * s;
        rx += length[k] * c; ry += length[k] * s;
    }
    const double det = gxx * gyy - gxy * gxy;
    if (!(det > 0)) fail(op, no_solution, "the conformal chart's boundary cannot be closed");
    const double lx = (gyy * rx - gxy * ry) / det, ly = (gxx * ry - gxy * rx) / det;
    std::vector<double> x(n, 0), y(n, 0);
    gp_XY at(0, 0);
    for (size_t k = 0; k < m; ++k) {
        const double c = std::cos(heading[k]), s = std::sin(heading[k]);
        x[chart.boundary[k]] = at.X();
        y[chart.boundary[k]] = at.Y();
        const double closed = length[k] * (1 - c * lx - s * ly);
        if (!(closed > 0)) fail(op, no_solution, "the conformal chart's boundary cannot be closed");
        at += gp_XY(c, s) * closed;
    }
    // The inside follows harmonically: a conformal map's coordinates are harmonic.
    const std::vector<double> none(n, 0);
    solve(laplacian, onBoundary, none, x);
    solve(laplacian, onBoundary, none, y);
    chart.flat.resize(n);
    for (int i = 0; i < n; ++i) chart.flat[i] = gp_XY(x[i], y[i]);

    // Centred on its area and turned to its principal axes.
    double area = 0, sx = 0, sy = 0;
    for (const auto &t : chart.triangles) {
        const double a = signedArea(chart.flat[t[0]], chart.flat[t[1]], chart.flat[t[2]]);
        if (!(a > 0)) fail(op, no_solution, "the conformal chart folds a triangle over");
        const gp_XY centre = (chart.flat[t[0]] + chart.flat[t[1]] + chart.flat[t[2]]) / 3;
        area += a; sx += a * centre.X(); sy += a * centre.Y();
    }
    const gp_XY centre(sx / area, sy / area);
    double cxx = 0, cxy = 0, cyy = 0;
    for (const auto &t : chart.triangles) {
        const double a = signedArea(chart.flat[t[0]], chart.flat[t[1]], chart.flat[t[2]]);
        for (int k = 0; k < 3; ++k) {
            const gp_XY d = chart.flat[t[k]] - centre;
            cxx += a * d.X() * d.X(); cxy += a * d.X() * d.Y(); cyy += a * d.Y() * d.Y();
        }
    }
    const double turn = std::atan2(2 * cxy, cxx - cyy) / 2, c = std::cos(turn), s = std::sin(turn);
    for (auto &p : chart.flat) {
        const gp_XY d = p - centre;
        p = gp_XY(c * d.X() + s * d.Y(), -s * d.X() + c * d.Y());
    }

    // The grid that finds a triangle by face parameters.
    gp_XY high = chart.source.front();
    chart.low = high;
    for (const auto &uv : chart.source) {
        chart.low.SetCoord(std::min(chart.low.X(), uv.X()), std::min(chart.low.Y(), uv.Y()));
        high.SetCoord(std::max(high.X(), uv.X()), std::max(high.Y(), uv.Y()));
    }
    const int side = std::max(1, static_cast<int>(std::sqrt(chart.triangles.size() / 2.0)));
    chart.columns = chart.rows = side;
    chart.cell = gp_XY((high.X() - chart.low.X()) / side, (high.Y() - chart.low.Y()) / side);
    if (!(chart.cell.X() > 0 && chart.cell.Y() > 0)) fail(op, no_solution, "the face to be flattened has no extent");
    chart.cells.assign(static_cast<size_t>(side) * side, {});
    const auto index = [&](double v, double low, double size) {
        return std::clamp(static_cast<int>(std::floor((v - low) / size)), 0, side - 1);
    };
    for (size_t k = 0; k < chart.triangles.size(); ++k) {
        const auto &t = chart.triangles[k];
        double u0 = HUGE_VAL, u1 = -HUGE_VAL, v0 = HUGE_VAL, v1 = -HUGE_VAL;
        for (int i : t) {
            u0 = std::min(u0, chart.source[i].X()); u1 = std::max(u1, chart.source[i].X());
            v0 = std::min(v0, chart.source[i].Y()); v1 = std::max(v1, chart.source[i].Y());
        }
        for (int a = index(u0, chart.low.X(), chart.cell.X()); a <= index(u1, chart.low.X(), chart.cell.X()); ++a)
            for (int b = index(v0, chart.low.Y(), chart.cell.Y()); b <= index(v1, chart.low.Y(), chart.cell.Y()); ++b)
                chart.cells[static_cast<size_t>(a) * side + b].push_back(static_cast<int>(k));
    }
    return chart;
}

gp_XY ConformalChart::map(const gp_XY &uv) const {
    const int a = std::clamp(static_cast<int>(std::floor((uv.X() - low.X()) / cell.X())), 0, columns - 1),
              b = std::clamp(static_cast<int>(std::floor((uv.Y() - low.Y()) / cell.Y())), 0, rows - 1);
    // The triangle nearest the place in the face's parameters: the one holding it, or for a
    // place just outside, the one whose edge it is nearest to.
    double best = HUGE_VAL;
    const std::array<int, 3> *found = nullptr;
    for (int ring = 0; ring <= std::max(columns, rows); ++ring) {
        for (int i = a - ring; i <= a + ring; ++i)
            for (int j = b - ring; j <= b + ring; ++j) {
                if (i < 0 || j < 0 || i >= columns || j >= rows || (std::abs(i - a) != ring && std::abs(j - b) != ring))
                    continue;
                for (int k : cells[static_cast<size_t>(i) * rows + j]) {
                    const auto &t = triangles[k];
                    double distance = 0;
                    for (int e = 0; e < 3; ++e) {
                        const gp_XY &p = source[t[e]], &q = source[t[(e + 1) % 3]];
                        if (signedArea(p, q, uv) >= 0) continue;
                        const gp_XY edge = q - p, to = uv - p;
                        distance = std::max(distance,
                            (to - edge * std::clamp(to.Dot(edge) / edge.SquareModulus(), 0.0, 1.0)).Modulus());
                    }
                    if (distance < best) {
                        best = distance;
                        found = &t;
                    }
                }
            }
        // Inside a triangle, or the nearest of those within a ring all round.
        if (found && (best == 0 || ring > 0)) break;
    }
    if (!found) fail(op, no_solution, "a place lies off the conformal chart");
    const auto &t = *found;
    const double whole = signedArea(source[t[0]], source[t[1]], source[t[2]]);
    return flat[t[0]] * (signedArea(uv, source[t[1]], source[t[2]]) / whole) +
           flat[t[1]] * (signedArea(source[t[0]], uv, source[t[2]]) / whole) +
           flat[t[2]] * (signedArea(source[t[0]], source[t[1]], uv) / whole);
}

} // namespace loft
