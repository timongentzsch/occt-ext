// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Phase.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_Band.hxx>
namespace BRepExt { using GeomAlgoExt::Band; }
#include <GeomAlgoExt_Failure.hxx>

#include <BSplCLib.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <math.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <algorithm>
#include <cmath>
#include <map>

#if !defined(__wasm_simd128__) || defined(OCCT_GEOMALGOEXT_SCALAR)
namespace BRepExt {
using GeomAlgoExt::v128_t;
using GeomAlgoExt::wasm_f64x2_splat;
using GeomAlgoExt::wasm_v128_load;
using GeomAlgoExt::wasm_v128_store;
using GeomAlgoExt::wasm_f64x2_add;
using GeomAlgoExt::wasm_f64x2_mul;
}
#endif

namespace BRepExt {
namespace {

/// What is added to the diagonal of the scaled matrix, whose entries are one, when it does not
/// factor: from the rounding of a sum of some hundred products of doubles (1e-13) in steps of
/// a hundred up to a millionth, beyond which the pull would show against the conditions.
constexpr double firstHold = 1e-13, holdStep = 100, lastHold = 1e-6;

/// The multipliers are given up when, over this many passes, the misses of the conditions have
/// fallen by less than this share: what is left then is what the net cannot follow.
constexpr int stallPasses = 12;

constexpr double stallShare = 0.02;

/// The most passes of the multipliers: with every one the miss falls or the stall above ends
/// them, and a net that is not there after some hundreds is not brought there by more.
constexpr int mostPasses = 400;

/// The multipliers have nothing left to work on when what drives them has fallen to this share
/// of what it was at first: the square of the rounding of a double.
constexpr double spent = 1e-20;

/// A condition over the unknowns it touches, scaled by the root of its weight; what the fixed
/// poles contribute is taken to its value.
struct Row {
    std::vector<int> index;
    std::vector<double> coefficient;
    double value, within;
    /// Which of the conditions it is.
    int from;
};

} // namespace

SplineBasis splineBasis(const char *op, const std::vector<double> &knots, int degree, double t) {
    const NCollection_Array1<double> flat(knots.front(), 1, static_cast<int>(knots.size()));
    math_Matrix values(1, splineDerivatives + 1, 1, degree + 1, 0.0);
    SplineBasis basis{};
    if (BSplCLib::EvalBsplineBasis(std::min(splineDerivatives, degree), degree + 1, flat, t, basis.first, values))
        fail(op, no_solution, "the patch's knots are degenerate");
    basis.first -= 1;
    for (int d = 0; d <= splineDerivatives; ++d)
        for (int k = 0; k <= degree; ++k)
            basis.value[d][k] = values(d + 1, k + 1);
    return basis;
}

std::vector<int> fitNet(const char *op, SplineNet &net, const std::vector<char> &fixed,
            const std::vector<FitCondition> &conditions, const FitEnergy &energy,
            const std::function<double(double, double)> &density, const gp_XYZ *along,
            const std::vector<FitJump> &jumps, bool leastSquares, FitEnergies *again) {
    PatchPhase phase("fit");
    if (leastSquares && !along)
        fail(op, no_solution, "the initial fair fit requires a height-graph direction");
    const int nu = net.countU(), nv = net.countV(), pu = net.degreeU, pv = net.degreeV;
    // Unknowns per pole: its three coordinates, or how far it stands along the one direction.
    const int dim = along ? 1 : 3;
    const int size = dim * nu * nv;
    // The energies and the jumps never mix the coordinates, and where no condition does either
    // each coordinate is a system of its own: they are then taken one after the other, each in
    // a band a third as wide, and not side by side pole by pole. Their entries and the order
    // of every sum are the same; only the zeros between them are not visited.
    const bool separate = dim == 3 && std::all_of(conditions.begin(), conditions.end(),
        [](const FitCondition &condition) {
            int axes = 0;
            for (int t = 0; t < condition.count; ++t)
                for (int c = 0; c < 3; ++c)
                    if (condition.terms[t].direction.Coord(c + 1) != 0) axes |= 1 << c;
            return (axes & (axes - 1)) == 0;
        });
    const bool reused = again && again->matrix && again->separate == separate && again->matrix->size == size;
    const auto unknown = [&](int p, int c) { return separate ? c * nu * nv + p : dim * p + c; };
    // A jump is read on both sides of a knot, and so over the poles of two spans.
    const int reach = jumps.empty() ? pu * nv + pv : 2 * pu * nv + 2 * pv;
    Band matrix(size, separate ? reach : dim * reach + dim - 1, !jumps.empty());
    matrix.threeBlocks = separate;
    std::vector<double> constant(size, 0);
    // What a vector counts for in unknown `c` of a pole, and a pole as it was given in it.
    const auto part = [&](const gp_XYZ &vector, int c) { return along ? vector.Dot(*along) : vector.Coord(c + 1); };

    // The energies, by Gauss points in every span: as many each way as integrate a product of
    // two derivatives exactly. The squared second and third derivatives by the parameters,
    // each mixed one as often as it occurs; the coordinates do not meet in them.
    const auto spans = [](const std::vector<double> &knots, int degree, std::vector<std::pair<double, double>> &out) {
        math_Vector points(1, degree + 1), weights(1, degree + 1);
        math::GaussPoints(degree + 1, points);
        math::GaussWeights(degree + 1, weights);
        for (size_t k = degree; k + degree + 1 < knots.size(); ++k)
            for (int g = 1; knots[k + 1] > knots[k] && g <= degree + 1; ++g)
                out.emplace_back((knots[k] + knots[k + 1]) / 2 + (knots[k + 1] - knots[k]) / 2 * points(g),
                                 (knots[k + 1] - knots[k]) / 2 * weights(g));
    };
    std::vector<std::pair<double, double>> us, vs;
    if (!reused) {
        spans(net.knotsU, pu, us);
        spans(net.knotsV, pv, vs);
    }
    std::vector<SplineBasis> bus, bvs;
    for (const auto &[u, weight] : us)
        bus.push_back(splineBasis(op, net.knotsU, pu, u));
    for (const auto &[v, weight] : vs)
        bvs.push_back(splineBasis(op, net.knotsV, pv, v));
    static const struct {
        int du, dv;
        double times;
    } parts[] = {{2, 0, 1}, {1, 1, 2}, {0, 2, 1}, {3, 0, 1}, {2, 1, 3}, {1, 2, 3}, {0, 3, 1}};
    const int local = (pu + 1) * (pv + 1);
    std::vector<double> element(static_cast<size_t>(local) * local), values(local);
    const auto addElement = [&](const SplineBasis &bu, const SplineBasis &bv) {
        auto pole = [&](int i) { return (bu.first + i / (pv + 1)) * nv + bv.first + i % (pv + 1); };
        for (int i = 0; i < local; ++i)
            for (int k = 0; k <= i; ++k) {
                const int p = pole(i), q = pole(k);
                const double e = element[static_cast<size_t>(i) * local + k];
                for (int c = 0; c < dim; ++c) {
                    // What is looked for is how far each pole moves from where it was
                    // given: what the given net bends goes to the other side.
                    // Separate coordinates have one energy: it is added to the first.
                    if (!fixed[p] && !fixed[q] && c == 0)
                        matrix.add(unknown(p, c), unknown(q, c), e);
                    const double givenP = part(net.poles[p], c), givenQ = part(net.poles[q], c);
                    if (!fixed[p])
                        constant[unknown(p, c)] -= e * givenQ;
                    if (!fixed[q] && i != k)
                        constant[unknown(q, c)] -= e * givenP;
                }
            }
    };
    // Independent SIMD lanes retain each entry's multiplication and addition order.
    const auto accumulate = [&](double scale) {
        for (int i = 0; i < local; ++i) {
            double *row = element.data() + static_cast<size_t>(i) * local;
            const v128_t weight = wasm_f64x2_splat(scale * values[i]);
            int k = 0;
            for (; k + 1 <= i; k += 2)
                wasm_v128_store(row + k, wasm_f64x2_add(wasm_v128_load(row + k),
                    wasm_f64x2_mul(weight, wasm_v128_load(values.data() + k))));
            if (k <= i) row[k] += scale * values[i] * values[k];
        }
    };
    PatchPhase assembly("assembly");
    const bool derivativeConstraints = std::any_of(conditions.begin(), conditions.end(),
        [](const FitCondition &condition) {
            for (int t = 0; t < condition.count; ++t)
                if (condition.terms[t].du + condition.terms[t].dv > 0) return true;
            return false;
        });
    // Free-curve fits can be ill-conditioned enough that quadrature summation order
    // changes the refinement decisions. Preserve their arithmetic; derivative-held
    // quads use the separable integration already used by corner-jump fits.
    // Fitted again from the same poles with other conditions: the energies as they were.
    if (reused) {
        matrix = *again->matrix;
        constant = again->constant;
    } else if (jumps.empty() && (!derivativeConstraints || density)) {
        for (size_t gu = 0; gu < us.size(); ++gu)
            for (size_t gv = 0; gv < vs.size(); ++gv) {
                const double here = us[gu].second * vs[gv].second * (density ? density(us[gu].first, vs[gv].first) : 1);
                const SplineBasis &bu = bus[gu], &bv = bvs[gv];
                std::fill(element.begin(), element.end(), 0);
                for (const auto &term : parts) {
                    for (int a = 0; a <= pu; ++a)
                        for (int b = 0; b <= pv; ++b)
                            values[a * (pv + 1) + b] = bu.value[term.du][a] * bv.value[term.dv][b];
                    const double scale =
                        here * term.times * (term.du + term.dv == 2 ? energy.bending : energy.variation);
                    accumulate(scale);
                }
                addElement(bu, bv);
            }
    } else {
        // Integrate the tensor factors once per span. Most cells have one density throughout;
        // only cells crossed by the trim need the full two-dimensional quadrature.
        const auto grams = [](const std::vector<SplineBasis> &basis,
                              const std::vector<std::pair<double, double>> &sites, int degree) {
            const int q = degree + 1;
            std::vector<std::vector<double>> out(sites.size() / q, std::vector<double>(4 * q * q, 0));
            for (size_t g = 0; g < sites.size(); ++g)
                for (int derivative = 0; derivative <= splineDerivatives; ++derivative)
                    for (int a = 0; a < q; ++a)
                        for (int b = 0; b < q; ++b)
                            out[g / q][(derivative * q + a) * q + b] +=
                                sites[g].second * basis[g].value[derivative][a] * basis[g].value[derivative][b];
            return out;
        };
        const auto gramU = grams(bus, us, pu), gramV = grams(bvs, vs, pv);
        const int qu = pu + 1, qv = pv + 1;
        std::vector<double> densities(qu * qv);
        for (size_t gu = 0; gu < us.size(); gu += qu)
            for (size_t gv = 0; gv < vs.size(); gv += qv) {
                const SplineBasis &bu = bus[gu], &bv = bvs[gv];
                bool uniform = true;
                for (int a = 0; a < qu; ++a)
                    for (int b = 0; b < qv; ++b) {
                        const double value = density ? density(us[gu + a].first, vs[gv + b].first) : 1;
                        densities[a * qv + b] = value;
                        uniform = uniform && value == densities[0];
                    }
                std::fill(element.begin(), element.end(), 0);
                if (uniform) {
                    const auto &u = gramU[gu / qu], &v = gramV[gv / qv];
                    for (int i = 0; i < local; ++i)
                        for (int k = 0; k <= i; ++k) {
                            double sum = 0;
                            for (const auto &term : parts)
                                sum += term.times * (term.du + term.dv == 2 ? energy.bending : energy.variation) *
                                    u[(term.du * qu + i / qv) * qu + k / qv] *
                                    v[(term.dv * qv + i % qv) * qv + k % qv];
                            element[static_cast<size_t>(i) * local + k] = densities[0] * sum;
                        }
                } else {
                    for (int a = 0; a < qu; ++a)
                        for (int b = 0; b < qv; ++b) {
                            const double here = us[gu + a].second * vs[gv + b].second * densities[a * qv + b];
                            for (const auto &term : parts) {
                                for (int i = 0; i < qu; ++i)
                                    for (int j = 0; j < qv; ++j)
                                        values[i * qv + j] = bus[gu + a].value[term.du][i] * bvs[gv + b].value[term.dv][j];
                                const double scale = here * term.times *
                                    (term.du + term.dv == 2 ? energy.bending : energy.variation);
                                accumulate(scale);
                            }
                        }
                }
                addElement(bu, bv);
            }
    }
    // Energy is scalar in every spatial coordinate, even when later conditions mix
    // them. Copy the completed sum, retaining its order and the original skyline.
    if (dim == 3 && !reused) {
        for (int p = 0; p < nu * nv; ++p) {
            const int row = unknown(p, 0), first = matrix.first[row];
            for (int c = 1; c < dim; ++c) {
                const int next = unknown(p, c), offset = next - row;
                matrix.first[next] = first + offset;
                for (int q = separate ? first : first / dim; q <= (separate ? row : p); ++q) {
                    const int column = separate ? q : unknown(q, 0);
                    matrix.row(next)[column + offset] = matrix.row(row)[column];
                }
            }
        }
    }
    for (int p = 0; !reused && p < nu * nv; ++p)
        if (fixed[p])
            for (int c = 0; c < dim; ++c) {
                matrix.add(unknown(p, c), unknown(p, c), 1);
                constant[unknown(p, c)] = 0;
            }
    if (again && !reused) {
        again->matrix = matrix;
        again->constant = constant;
        again->separate = separate;
    }

    std::vector<Row> rows, watched;
    rows.reserve(conditions.size());
    for (const auto &condition : conditions) {
        const auto bu = splineBasis(op, net.knotsU, pu, condition.u), bv = splineBasis(op, net.knotsV, pv, condition.v);
        Row row;
        row.from = static_cast<int>(&condition - conditions.data());
        const double scale = std::sqrt(condition.weight);
        row.value = condition.value * scale;
        row.within = condition.within * scale;
        for (int a = 0; a <= pu; ++a)
            for (int b = 0; b <= pv; ++b) {
                const int p = (bu.first + a) * nv + bv.first + b;
                gp_XYZ direction(0, 0, 0);
                for (int t = 0; t < condition.count; ++t)
                    direction += condition.terms[t].direction *
                                 (bu.value[condition.terms[t].du][a] * bv.value[condition.terms[t].dv][b]);
                row.value -= direction.Dot(net.poles[p]) * scale;
                for (int c = 0; c < dim && !fixed[p]; ++c)
                    if (const double coefficient = part(direction, c) * scale; coefficient != 0) {
                        row.index.push_back(unknown(p, c));
                        row.coefficient.push_back(coefficient);
                    }
            }
        // Nothing to be done about it by what is free: a point's place across the direction.
        if (row.index.empty())
            continue;
        if (condition.watch) {
            watched.push_back(std::move(row));
            continue;
        }
        for (size_t a = 0; a < row.index.size(); ++a)
            for (size_t b = 0; b <= a; ++b)
                matrix.add(row.index[a], row.index[b], row.coefficient[a] * row.coefficient[b]);
        if (condition.soft) {
            for (size_t a = 0; a < row.index.size(); ++a)
                constant[row.index[a]] += row.value * row.coefficient[a];
            continue;
        }
        rows.push_back(std::move(row));
    }
    // Reversing the knots reads the exact left limit at a repeated knot. Reading a
    // nearby parameter instead introduces tiny coefficients on otherwise absent poles,
    // widening the Cholesky skyline and every multiplier solve for no useful constraint.
    const auto leftBasis = [&](const std::vector<double> &knots, int degree, double t) {
        std::vector<double> reversed(knots.rbegin(), knots.rend());
        for (double &knot : reversed) knot = -knot;
        SplineBasis basis = splineBasis(op, reversed, degree, -t);
        basis.first = static_cast<int>(knots.size()) - 2 * degree - 2 - basis.first;
        for (int d = 0; d <= splineDerivatives; ++d) {
            std::reverse(basis.value[d], basis.value[d] + degree + 1);
            if (d % 2)
                for (int k = 0; k <= degree; ++k) basis.value[d][k] = -basis.value[d][k];
        }
        return basis;
    };
    for (const FitJump &jump : jumps) {
        std::map<int, double> over;
        const auto &acrossKnots = jump.acrossU ? net.knotsU : net.knotsV;
        const int acrossDegree = jump.acrossU ? pu : pv;
        const auto along = splineBasis(op, jump.acrossU ? net.knotsV : net.knotsU,
                                       jump.acrossU ? pv : pu, jump.along);
        for (const double side : {1.0, -1.0}) {
            const auto across = side > 0 ? splineBasis(op, acrossKnots, acrossDegree, jump.at)
                                         : leftBasis(acrossKnots, acrossDegree, jump.at);
            const auto &bu = jump.acrossU ? across : along, &bv = jump.acrossU ? along : across;
            for (int a = 0; a <= pu; ++a)
                for (int b = 0; b <= pv; ++b)
                    over[(bu.first + a) * nv + bv.first + b] +=
                        side * bu.value[jump.acrossU ? 2 : 0][a] * bv.value[jump.acrossU ? 0 : 2][b];
        }
        const double scale = std::sqrt(jump.weight);
        for (int c = 0; c < 3; ++c) {
            gp_XYZ axis(0, 0, 0);
            axis.SetCoord(c + 1, 1);
            Row row;
            row.value = 0;
            for (const auto &[p, weight] : over) {
                row.value -= weight * axis.Dot(net.poles[p]) * scale;
                for (int k = 0; k < dim && !fixed[p]; ++k)
                    if (const double coefficient = part(axis, k) * weight * scale; coefficient != 0) {
                        // One unknown once: along one direction the three coordinates share it.
                        const int index = unknown(p, k);
                        const auto found = std::find(row.index.begin(), row.index.end(), index);
                        if (found == row.index.end()) {
                            row.index.push_back(index);
                            row.coefficient.push_back(coefficient);
                        } else {
                            row.coefficient[found - row.index.begin()] += coefficient;
                        }
                    }
            }
            for (size_t a = 0; a < row.index.size(); ++a) {
                for (size_t b = 0; b <= a; ++b)
                    matrix.add(row.index[a], row.index[b], row.coefficient[a] * row.coefficient[b]);
                constant[row.index[a]] += row.value * row.coefficient[a];
            }
        }
    }
    // The weights of the conditions and of the energies are many orders apart, and so are the
    // rows of the matrix. Each unknown is taken in the unit that makes its own diagonal entry
    // one (the scaling of Jacobi), so that what the factorisation sees is how the unknowns hang
    // together and not how they were weighed.
    std::vector<double> unit(size, 1);
    for (int i = 0; i < size; ++i)
        if (const double diagonal = matrix.row(i)[i]; diagonal > 0 && std::isfinite(diagonal))
            unit[i] = 1 / std::sqrt(diagonal);
    for (int i = 0; i < size; ++i) {
        double *const ri = matrix.row(i);
        for (int k = std::max(0, i - matrix.half); k <= i; ++k)
            ri[k] *= unit[i] * unit[k];
        constant[i] *= unit[i];
    }
    for (auto *list : {&rows, &watched})
        for (Row &row : *list)
            for (size_t k = 0; k < row.index.size(); ++k)
                row.coefficient[k] *= unit[row.index[k]];
    assembly.finish();
    // Poles that no condition reaches and the energy holds only by rounding (far outside a
    // loop, at a thousandth of a weight) can still leave the matrix short of definite. They
    // are then held, as weakly as it takes, to where they were given: a pull that is nothing
    // against any condition and decides only what nothing else does.
    {
        PatchPhase phase("factorization");
        const Band assembled = matrix;
        double hold = 0;
        while (!matrix.factor()) {
            hold = hold == 0 ? firstHold : hold * holdStep;
            if (hold > lastHold)
                fail(op, no_solution,
                     "the equations of the fit cannot be solved: the net is degenerate, or a "
                     "condition is not a number");
            matrix = assembled;
            for (int i = 0; i < size; ++i)
                matrix.row(i)[i] += hold;
        }
    }

    // The conditions are weighted, not imposed, so the solution of this system leaves each a
    // little short by what the energy pulls against it. A shift of each condition's value (a
    // multiplier) takes that away, down to what the net cannot follow at all. The shifts solve
    // a symmetric system of their own, whose product with a vector is one solve with the factor
    // above; it is singular where conditions say the same thing twice, which the method of
    // conjugate residuals does not mind.
    const size_t count = rows.size();
    // The unknowns for given shifts, and what the conditions then miss by.
    auto poles = [&](const std::vector<double> &shift, bool homogeneous) {
        std::vector<double> x = homogeneous ? std::vector<double>(size, 0) : constant;
        for (size_t r = 0; r < count; ++r) {
            const double pull = homogeneous ? shift[r] : rows[r].value - shift[r];
            for (size_t k = 0; k < rows[r].index.size(); ++k)
                x[rows[r].index[k]] += pull * rows[r].coefficient[k];
        }
        { PatchPhase phase("linear-solves"); matrix.solve(x); }
        return x;
    };
    auto apply = [&](const std::vector<double> &x, bool homogeneous) {
        std::vector<double> out(count);
        for (size_t r = 0; r < count; ++r) {
            double sum = homogeneous ? 0 : -rows[r].value;
            for (size_t k = 0; k < rows[r].index.size(); ++k)
                sum += rows[r].coefficient[k] * x[rows[r].index[k]];
            out[r] = sum;
        }
        return out;
    };
    auto dot = [](const std::vector<double> &a, const std::vector<double> &b) {
        double sum = 0;
        for (size_t i = 0; i < a.size(); ++i)
            sum += a[i] * b[i];
        return sum;
    };
    std::vector<double> shift(count, 0);
    // The unknowns are kept in step with the shifts, for the conditions that are only read.
    std::vector<double> x = poles(shift, false);
    const auto read = [&](const Row &row) {
        double sum = -row.value;
        for (size_t k = 0; k < row.index.size(); ++k)
            sum += row.coefficient[k] * x[row.index[k]];
        return std::abs(sum) <= row.within;
    };
    if (!leastSquares) {
        std::vector<double> residual = apply(x, false);
        std::vector<double> moved = poles(residual, true);
        std::vector<double> image = apply(moved, true), direction = residual, step = image;
        double energy0 = dot(residual, image);
        const double start = energy0;
        double left = dot(residual, residual), leftBefore = left;
        // Every condition as near as it was asked to be: nearer is not better, and costs the
        // surface some of its fairness.
        const auto met = [&] {
            for (size_t r = 0; r < count; ++r)
                if (std::abs(residual[r]) > rows[r].within)
                    return false;
            return std::all_of(watched.begin(), watched.end(), read);
        };
        for (int pass = 0; pass < mostPasses && energy0 > spent * start && !met(); ++pass) {
            const double alpha = energy0 / dot(step, step);
            for (size_t r = 0; r < count; ++r) {
                shift[r] += alpha * direction[r];
                residual[r] -= alpha * step[r];
            }
            for (int i = 0; i < size; ++i)
                x[i] -= alpha * moved[i];
            left = dot(residual, residual);
            if (pass % stallPasses == stallPasses - 1) {
                if (left > (1 - stallShare) * leftBefore)
                    break;
                leftBefore = left;
            }
            const std::vector<double> next = poles(residual, true);
            image = apply(next, true);
            const double energy1 = dot(residual, image), beta = energy1 / energy0;
            for (size_t r = 0; r < count; ++r) {
                direction[r] = residual[r] + beta * direction[r];
                step[r] = image[r] + beta * step[r];
            }
            for (int i = 0; i < size; ++i)
                moved[i] = next[i] + beta * moved[i];
            energy0 = energy1;
        }
        x = poles(shift, false);
    }
    if (leastSquares && !std::all_of(x.begin(), x.end(), [](double value) { return std::isfinite(value); }))
        fail(op, no_solution, "the initial fair fit is not finite");
    for (int p = 0; p < nu * nv; ++p) {
        if (fixed[p])
            continue;
        if (along)
            net.poles[p] += *along * (x[p] * unit[p]);
        else
            net.poles[p] += gp_XYZ(x[unknown(p, 0)] * unit[unknown(p, 0)], x[unknown(p, 1)] * unit[unknown(p, 1)],
                                   x[unknown(p, 2)] * unit[unknown(p, 2)]);
    }
    std::vector<int> missed;
    for (const Row &row : watched)
        if (!read(row))
            missed.push_back(row.from);
    return missed;
}

occ::handle<Geom_BSplineSurface> surfaceOf(const SplineNet &net) {
    auto distinct = [](const std::vector<double> &flat, std::vector<double> &knots, std::vector<int> &counts) {
        for (double knot : flat)
            if (!knots.empty() && knot == knots.back()) {
                ++counts.back();
            } else {
                knots.push_back(knot);
                counts.push_back(1);
            }
    };
    std::vector<double> ku, kv;
    std::vector<int> mu, mv;
    distinct(net.knotsU, ku, mu);
    distinct(net.knotsV, kv, mv);
    const int nu = net.countU(), nv = net.countV();
    NCollection_Array2<gp_Pnt> poles(1, nu, 1, nv);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
            poles.SetValue(i + 1, j + 1, gp_Pnt(net.poles[static_cast<size_t>(i) * nv + j]));
    return new Geom_BSplineSurface(poles, NCollection_Array1<double>(ku.front(), 1, static_cast<int>(ku.size())),
                                   NCollection_Array1<double>(kv.front(), 1, static_cast<int>(kv.size())),
                                   NCollection_Array1<int>(mu.front(), 1, static_cast<int>(mu.size())),
                                   NCollection_Array1<int>(mv.front(), 1, static_cast<int>(mv.size())),
                                   net.degreeU, net.degreeV);
}

} // namespace BRepExt
