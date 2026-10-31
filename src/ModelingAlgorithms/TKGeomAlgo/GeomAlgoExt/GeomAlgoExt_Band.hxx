// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once

#include <GeomAlgoExt_Pair.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

namespace GeomAlgoExt {

/// The lower half of a symmetric matrix whose entries lie within `half` of the diagonal.
struct Band {
    int size, half;
    std::vector<double> data;
    std::vector<int> first;
    bool paired;
    /// Three equal blocks along the diagonal that nothing joins, each block's rows beginning
    /// where the first's do: see sideBySide.
    bool threeBlocks = false;
    Band(int size, int half, bool paired)
        : size(size), half(half), data(static_cast<size_t>(size) * (half + 1), 0), first(size), paired(paired) {
        for (int i = 0; i < size; ++i) first[i] = i;
    }

    /// Whether the matrix is three equal blocks along the diagonal as threeBlocks says: then
    /// factor and solve take the same row of the three side by side, so that their sums, each
    /// formed in its own order as before, overlap.
    bool sideBySide() const {
        if (size % 3 || paired) return false;
        const int count = size / 3;
        for (int c = 1; c < 3; ++c)
            for (int i = 0; i < count; ++i)
                if (first[c * count + i] != first[i] + c * count) return false;
        return true;
    }
    double *row(int r) { return data.data() + static_cast<size_t>(r) * (half + 1) + half - r; }
    void add(int r, int c, double value) {
        if (value == 0) return;
        if (r < c) std::swap(r, c);
        row(r)[c] += value;
        first[r] = std::min(first[r], c);
    }

    // The supported fits with jump penalties spend most of their time in these products.
    // Pair their independent lanes; ordinary fits retain the original summation order.
    double subtract(double sum, const double *a, const double *b, int count) const {
        int k = 0;
        if (paired) {
            v128_t products = wasm_f64x2_splat(0), next = wasm_f64x2_splat(0);
            for (; k + 3 < count; k += 4) {
                products = wasm_f64x2_add(products,
                    wasm_f64x2_mul(wasm_v128_load(a + k), wasm_v128_load(b + k)));
                next = wasm_f64x2_add(next,
                    wasm_f64x2_mul(wasm_v128_load(a + k + 2), wasm_v128_load(b + k + 2)));
            }
            products = wasm_f64x2_add(products, next);
            if (k + 1 < count) {
                products = wasm_f64x2_add(products,
                    wasm_f64x2_mul(wasm_v128_load(a + k), wasm_v128_load(b + k)));
                k += 2;
            }
            sum -= wasm_f64x2_extract_lane(products, 0);
            sum -= wasm_f64x2_extract_lane(products, 1);
        }
        for (; k < count; ++k) sum -= a[k] * b[k];
        return sum;
    }

    /// Replace the matrix by its Cholesky factor; false when it is not positive definite.
    bool factor() {
        if (threeBlocks && sideBySide()) {
            const int count = size / 3;
            const size_t stride = static_cast<size_t>(count) * (half + 1);
            for (int i = 0; i < count; ++i) {
                double *const ri = row(i);
                const int from = first[i];
                for (int j = from; j <= i; ++j) {
                    const double *const rj = row(j);
                    const int begin = std::max(from, first[j]);
                    double sum[3] = {ri[j], ri[j + stride], ri[j + 2 * stride]};
                    for (int k = begin; k < j; ++k) {
                        sum[0] -= ri[k] * rj[k];
                        sum[1] -= ri[k + stride] * rj[k + stride];
                        sum[2] -= ri[k + 2 * stride] * rj[k + 2 * stride];
                    }
                    for (int c = 0; c < 3; ++c) {
                        if (j < i)
                            ri[j + c * stride] = sum[c] / rj[j + c * stride];
                        else if (sum[c] <= 0 || !std::isfinite(sum[c]))
                            return false;
                        else
                            ri[i + c * stride] = std::sqrt(sum[c]);
                    }
                }
            }
            return true;
        }
        for (int i = 0; i < size; ++i) {
            double *const ri = row(i);
            const int from = first[i];
            for (int j = from; j <= i; ++j) {
                const double *const rj = row(j);
                const int begin = std::max(from, first[j]);
                const double sum = subtract(ri[j], ri + begin, rj + begin, j - begin);
                if (j < i)
                    ri[j] = sum / rj[j];
                else if (sum <= 0 || !std::isfinite(sum))
                    return false;
                else
                    ri[i] = std::sqrt(sum);
            }
        }
        return true;
    }

    /// Solve in place with the factor.
    void solve(std::vector<double> &x) {
        if (threeBlocks && sideBySide()) {
            const int count = size / 3;
            const size_t stride = static_cast<size_t>(count) * (half + 1);
            double *const y[3] = {x.data(), x.data() + count, x.data() + 2 * count};
            for (int i = 0; i < count; ++i) {
                const double *const ri = row(i);
                double sum[3] = {y[0][i], y[1][i], y[2][i]};
                for (int k = first[i]; k < i; ++k) {
                    sum[0] -= ri[k] * y[0][k];
                    sum[1] -= ri[k + stride] * y[1][k];
                    sum[2] -= ri[k + 2 * stride] * y[2][k];
                }
                for (int c = 0; c < 3; ++c) y[c][i] = sum[c] / ri[i + c * stride];
            }
        } else
        for (int i = 0; i < size; ++i) {
            const double *const ri = row(i);
            const double sum = subtract(x[i], ri + first[i], x.data() + first[i], i - first[i]);
            x[i] = sum / ri[i];
        }
        for (int i = size - 1; i >= 0; --i) {
            x[i] /= row(i)[i];
            const double *const ri = row(i);
            int k = first[i];
            // Each entry takes the same steps as alone, two at a time: nothing is summed here.
            {
                const v128_t here = wasm_f64x2_splat(x[i]);
                const auto move = [&](int at) {
                    wasm_v128_store(x.data() + at, wasm_f64x2_sub(wasm_v128_load(x.data() + at),
                        wasm_f64x2_mul(wasm_v128_load(ri + at), here)));
                };
                for (; k + 3 < i; k += 4) { move(k); move(k + 2); }
                if (k + 1 < i) { move(k); k += 2; }
            }
            for (; k < i; ++k) x[k] -= ri[k] * x[i];
        }
    }
};

} // namespace GeomAlgoExt
