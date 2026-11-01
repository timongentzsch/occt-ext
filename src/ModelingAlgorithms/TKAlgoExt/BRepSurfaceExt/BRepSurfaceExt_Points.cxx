// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BSplCLib.hxx>
#include <algorithm>
#include <cmath>
#include <math_Matrix.hxx>

namespace BRepExt {
std::vector<double> surfacePoles(OcctKernel &k, uint32_t id) {
    const char *op = "surfacePoles";
    auto face = surfacing::faceOf(k, op, id);
    auto surface = surfacing::editable(op, face, false);
    std::vector<double> out{double(surface->NbUPoles()), double(surface->NbVPoles())};
    for (int i = 1; i <= surface->NbUPoles(); ++i)
        for (int j = 1; j <= surface->NbVPoles(); ++j) {
            auto p = surface->Pole(i, j);
            out.insert(out.end(), {p.X(), p.Y(), p.Z(), surface->Weight(i, j)});
        }
    const auto net = surfacing::netOf(surface);
    auto greville = [](const std::vector<double> &knots, int index, int degree) {
        double sum = 0;
        for (int j = 1; j <= degree; ++j) sum += knots[index + j];
        return sum / degree;
    };
    // The surface normal belongs to the pole's Greville parameters, not its hull triangle.
    for (int i = 0; i < surface->NbUPoles(); ++i)
        for (int j = 0; j < surface->NbVPoles(); ++j) {
            gp_Pnt p;
            gp_Vec du, dv;
            surface->D1(greville(net.knotsU, i, net.degreeU),
                        greville(net.knotsV, j, net.degreeV), p, du, dv);
            auto normal = du.Crossed(dv);
            if (normal.SquareMagnitude() > 1e-24) normal.Normalize();
            else normal = gp_Vec(0, 0, 0);
            out.insert(out.end(), {normal.X(), normal.Y(), normal.Z()});
        }
    return out;
}
namespace {
double dot(const std::vector<double> &a, const std::vector<double> &b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); ++i)
        sum += a[i] * b[i];
    return sum;
}
// Orthogonal projection onto the affine constraint space. The orthonormal rows span its
// complement: remove precisely that part of the requested pole offsets. Reorthogonalise to
// handle the redundant endpoint and tangent equations without a singular normal system.
void project(const char *op, SplineNet &net, const SplineNet &desired,
             const std::vector<FitCondition> &conditions) {
    const int count = 3 * net.poles.size();
    struct Row {
        std::vector<double> values;
        double rhs, within;
    };
    std::vector<Row> rows;
    std::vector<double> given;
    for (const auto &p : net.poles)
        given.insert(given.end(), {p.X(), p.Y(), p.Z()});
    for (const auto &condition : conditions) {
        math_Matrix bu(1, 3, 1, net.degreeU + 1, 0.0), bv(1, 3, 1, net.degreeV + 1, 0.0);
        int iu, iv;
        NCollection_Array1<double> ku(net.knotsU.front(), 1, net.knotsU.size()),
            kv(net.knotsV.front(), 1, net.knotsV.size());
        BSplCLib::EvalBsplineBasis(std::min(2, net.degreeU), net.degreeU + 1, ku, condition.u, iu, bu);
        BSplCLib::EvalBsplineBasis(std::min(2, net.degreeV), net.degreeV + 1, kv, condition.v, iv, bv);
        std::vector<double> row(count, 0);
        for (int i = 0; i <= net.degreeU; ++i)
            for (int j = 0; j <= net.degreeV; ++j) {
                auto direction = gp_XYZ(0, 0, 0);
                for (int t = 0; t < condition.count; ++t)
                    direction += condition.terms[t].direction * bu(condition.terms[t].du + 1, i + 1) *
                                 bv(condition.terms[t].dv + 1, j + 1);
                const int pole = (iu - 1 + i) * net.countV() + iv - 1 + j;
                for (int c = 0; c < 3; ++c)
                    row[3 * pole + c] = direction.Coord(c + 1);
            }
        const double norm = std::sqrt(dot(row, row));
        if (norm < 1e-14)
            continue;
        const double value = (condition.value - dot(row, given)) / norm;
        for (auto &c : row)
            c /= norm;
        rows.push_back({std::move(row), value, std::max(1e-9, condition.within / norm)});
    }
    std::vector<std::vector<double>> q;
    std::vector<double> rhs;
    // Pivot on the largest remaining row. Taking consecutive edge samples first makes the
    // first few rows almost parallel, which destroys rank decisions in streaming Gram-Schmidt.
    while (!rows.empty()) {
        size_t pivot = 0;
        double norm2 = 0;
        for (size_t i = 0; i < rows.size(); ++i) {
            double n = dot(rows[i].values, rows[i].values);
            if (n > norm2) {
                norm2 = n;
                pivot = i;
            }
        }
        if (norm2 < 1e-20) {
            for (const auto &row : rows)
                if (std::abs(row.rhs) > row.within)
                    fail(op, no_solution,
                         "neighbour constraints conflict; no pole displacement preserves them");
            break;
        }
        auto selected = std::move(rows[pivot]);
        rows.erase(rows.begin() + pivot);
        const double norm = std::sqrt(norm2);
        for (auto &c : selected.values)
            c /= norm;
        selected.rhs /= norm;
        // Reorthogonalise the pivot, and then each unused row, to machine precision.
        for (size_t i = 0; i < q.size(); ++i) {
            const double part = dot(selected.values, q[i]);
            for (int c = 0; c < count; ++c)
                selected.values[c] -= part * q[i][c];
            selected.rhs -= part * rhs[i];
        }
        const double unit = std::sqrt(dot(selected.values, selected.values));
        for (auto &c : selected.values)
            c /= unit;
        selected.rhs /= unit;
        for (auto &row : rows)
            for (int pass = 0; pass < 2; ++pass) {
                const double part = dot(row.values, selected.values);
                for (int c = 0; c < count; ++c)
                    row.values[c] -= part * selected.values[c];
                row.rhs -= part * selected.rhs;
            }
        q.push_back(std::move(selected.values));
        rhs.push_back(selected.rhs);
    }
    std::vector<double> x;
    for (size_t i = 0; i < desired.poles.size(); ++i) {
        auto d = desired.poles[i] - net.poles[i];
        x.insert(x.end(), {d.X(), d.Y(), d.Z()});
    }
    for (size_t i = 0; i < q.size(); ++i) {
        const double move = rhs[i] - dot(q[i], x);
        for (int c = 0; c < count; ++c)
            x[c] += move * q[i][c];
    }
    for (size_t i = 0; i < net.poles.size(); ++i)
        net.poles[i] += gp_XYZ(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
}
} // namespace
EvolutionData moveSurfacePoles(OcctKernel &k, uint32_t faceId, const std::vector<int> &indices,
                               const std::vector<double> &offsets, const std::vector<uint32_t> &edges,
                               const std::vector<uint32_t> &references, const std::vector<uint32_t> &supports,
                               const std::vector<int> &orders, double tolerance) {
    const char *op = "moveSurfacePoles";
    if (indices.empty() || indices.size() % 2 || offsets.size() != indices.size() / 2 * 3 ||
        edges.size() != references.size() || edges.size() != supports.size() ||
        edges.size() != orders.size() || edges.size() > 4 || !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input,
             "pole index pairs, one finite offset per pole, up to four neighbour constraints and a positive "
             "tolerance are required");
    const auto face = surfacing::faceOf(k, op, faceId);
    auto surface = surfacing::editable(op, face);
    auto net = surfacing::netOf(surface);
    if (net.poles.size() > 400)
        fail(op, unsupported_geometry, "continuity projection is limited to 400 poles");
    auto desired = net;
    std::vector<char> seen(net.poles.size(), 0);
    for (size_t i = 0; i < indices.size() / 2; ++i) {
        const int u = indices[2 * i], v = indices[2 * i + 1];
        if (u < 0 || v < 0 || u >= net.countU() || v >= net.countV())
            fail(op, invalid_input, "pole indices are zero-based and must lie inside the control net");
        if (seen[u * net.countV() + v])
            fail(op, invalid_input, "a pole is listed more than once");
        seen[u * net.countV() + v] = 1;
        if (!std::isfinite(offsets[3 * i]) || !std::isfinite(offsets[3 * i + 1]) ||
            !std::isfinite(offsets[3 * i + 2]))
            fail(op, invalid_input, "pole offsets must be finite");
        desired.poles[u * net.countV() + v] += gp_XYZ(offsets[3 * i], offsets[3 * i + 1], offsets[3 * i + 2]);
    }
    std::vector<surfacing::Constraint> targets;
    for (size_t i = 0; i < edges.size(); ++i) {
        auto own = surfacing::boundaryOf(k, op, edges[i], faceId),
             ref = surfacing::boundaryOf(k, op, references[i], supports[i]);
        if (face.IsSame(ref.face))
            fail(op, invalid_input, "a neighbour must be a different face");
        targets.push_back(surfacing::constraint(op, face, own.edge, ref, orders[i]));
    }
    // Keeping continuity presupposes it exists. Refuse before projecting an already broken joint.
    surfacing::measure(op, face, targets, tolerance);
    const bool g2 = std::any_of(orders.begin(), orders.end(), [](int order) { return order == 2; });
    for (int pass = 0; pass < (g2 ? 8 : 1); ++pass) {
        project(op, net, desired, surfacing::conditions(op, targets, surface, tolerance));
        surface = surfaceOf(net);
    }
    const auto made = surfacing::makeFace(op, surface, tolerance, face.Orientation());
    surfacing::measure(op, made, targets, tolerance);
    return stored(k, made);
}
} // namespace BRepExt
