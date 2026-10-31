// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#include <GeomAlgoExt_BlendFit.hxx>
#include <GeomAlgoExt_Execution.hxx>
#include <GeomAlgoExt_Hermite.hxx>
#include <BSplCLib.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <math_Matrix.hxx>
#include <algorithm>
#include <cmath>
namespace GeomAlgoExt {
static std::vector<double> flatOf(const std::map<double,int>& knots) {
 std::vector<double> flat;
 for (const auto& [u,m] : knots) flat.insert(flat.end(),m,u);
 return flat;
}
BlendSectionFitter::BlendSectionFitter(const BlendFitInput& input,std::map<double,int> knots)
 : myInput(input), myKnots(std::move(knots)) {}
bool BlendSectionFitter::Advance(int refinement) {
 const auto& input=myInput;
 auto& worst=myWorst;auto& bendRounds=myBendRounds;
 const int kDegree=input.degree, S=input.sectionCount, width=4*S+4;
 const double tolerance=input.tolerance;
 auto& knots=myKnots;
 const auto& row=input.row;
 const auto& terminalRow=input.terminalRow;
 const bool keepFirst=input.keepFirst,keepLast=input.keepLast,contactPilot=input.contactPilot;
 struct { int degree; const std::vector<double>& knots; const std::vector<int>& multiplicities; }
 reference{input.sectionDegree,input.sectionKnots,input.sectionMultiplicities};
  // Breaks found while refining, on the knot they fall on where they fall on one.
  const auto raise = [&](const std::map<double, int> &breaks) {
    bool more = false;
    for (const auto &[u, m] : breaks) {
      const auto near = knots.lower_bound(u - input.knotMerge);
      int &slot = near != knots.end() && near->first <= u + input.knotMerge ? near->second : knots[u];
      if (slot < m) {
        slot = m;
        more = true;
      }
    }
    return more;
  };

    const auto flat = flatOf(knots);
    const int count = int(flat.size()) - kDegree - 1;
    if (count > 4000) { myStopped=true; return false; }
    NCollection_Array1<double> flatKnots(1, int(flat.size()));
    for (size_t k = 0; k < flat.size(); ++k)
      flatKnots(int(k) + 1) = flat[k];
    std::vector<double> sites(count);
    for (int i = 0; i < count; ++i) {
      double g = 0;
      for (int j = 1; j <= kDegree; ++j)
        g += flat[i + j];
      sites[i] = g / kDegree;
    }
    std::vector<GeomAlgoExt::HermiteSite> jetSites;
    if (keepFirst || keepLast) {
      jetSites =
          GeomAlgoExt::hermiteSites(kDegree, flat, keepFirst ? 2 : 0, keepLast ? 2 : 0);
      if (jetSites.empty())
        fail("blend-fit", no_solution, "invalid terminal Hermite sites");
    }
    // Interpolation at the sites, all rows at once: the collocation matrix is
    // banded and, at Greville value sites (or terminal Hermite jets),
    // OpenCascade's banded elimination without pivoting solves it in linear
    // time.
    const auto basisAt = [&](double u, int &firstIndex, math_Matrix &values) {
      return BSplCLib::EvalBsplineBasis(0, kDegree + 1, flatKnots, u,
                                        firstIndex, values);
    };
    NCollection_Array1<double> parameters(1, count);
    NCollection_Array1<int> orders(1, count);
    std::vector<double> poles(size_t(count) * width);
    for (int i = 0; i < count; ++i) {
      const GeomAlgoExt::HermiteSite site =
          jetSites.empty() ? GeomAlgoExt::HermiteSite{sites[i], 0} : jetSites[i];
      parameters(i + 1) = site.parameter;
      orders(i + 1) = site.order;
      const bool terminal = (keepFirst && site.parameter == 0) ||
                            (keepLast && site.parameter == 1);
      const auto r = terminal ? terminalRow(site.parameter, site.order)
                              : row(site.parameter);
      std::copy(r.begin(), r.end(), poles.begin() + size_t(i) * width);
    }
    int singular = 0;
    BSplCLib::Interpolate(kDegree, flatKnots, parameters, orders, width,
                          poles[0], singular);
    if (singular != 0)
      fail("blend-fit", no_solution, "the blend's interpolation is singular");
    std::vector<std::vector<double>> coefficients(width,
                                                  std::vector<double>(count));
    for (int c = 0; c < width; ++c)
      for (int i = 0; i < count; ++i)
        coefficients[c][i] = poles[size_t(i) * width + c];
    const auto evaluate = [&](double u) {
      math_Matrix values(1, 1, 1, kDegree + 1);
      int firstIndex = 0;
      basisAt(u, firstIndex, values);
      std::vector<double> r(width, 0.0);
      for (int c = 0; c < width; ++c)
        for (int k = 1; k <= kDegree + 1; ++k) {
          const int i = firstIndex + k - 2;
          if (i >= 0 && i < count)
            r[c] += values(1, k) * coefficients[c][i];
        }
      return r;
    };
    // A G2 blend's normal curvature along each contact, against its support's
    // along the feet: the fit is G2 to its supports only so far as its contact
    // rows are, which the positions alone do not bound (beside a C1 knot where
    // an edge's parameter kinks, they bend apart by 3e-4 / r within the length
    // tolerance).
    const auto alongBend = [&](double u, int side) {
      math_Matrix values(1, 3, 1, kDegree + 1);
      int firstIndex = 0;
      BSplCLib::EvalBsplineBasis(2, kDegree + 1, flatKnots, u, firstIndex,
                                 values);
      const auto row = [&](int c, int d) {
        double sum = 0;
        for (int k = 1; k <= kDegree + 1; ++k) {
          const int i = firstIndex + k - 2;
          if (i >= 0 && i < count)
            sum += values(d + 1, k) * coefficients[c][i];
        }
        return sum;
      };
      const int j = side == 0 ? 0 : S - 1;
      gp_XYZ X[3];
      double w[3];
      for (int d = 0; d < 3; ++d) {
        X[d] = gp_XYZ(row(4 * j, d), row(4 * j + 1, d), row(4 * j + 2, d));
        w[d] = row(4 * j + 3, d);
      }
      const gp_XYZ P = X[0] / w[0], P1 = (X[1] - P * w[1]) / w[0],
                   P2 = (X[2] - P1 * (2 * w[1]) - P * w[2]) / w[0];
      const double U = row(4 * S + 2 * side, 0),
                   V = row(4 * S + 2 * side + 1, 0);
      const double U1 = row(4 * S + 2 * side, 1),
                   V1 = row(4 * S + 2 * side + 1, 1);
      BlendContactJet jet;
      if (!input.contactJet(side, U, V, jet))
        return 0.0;
      const gp_XYZ n = jet.n.XYZ(), c1 = (jet.du * U1 + jet.dv * V1).XYZ();
      const gp_XYZ c2 =
          (jet.duu * (U1 * U1) + jet.duv * (2 * U1 * V1) + jet.dvv * (V1 * V1))
              .XYZ();
      if (!(P1.SquareModulus() > 0 && c1.SquareModulus() > 0))
        return 0.0;
      return std::abs(n.Dot(P2) / P1.SquareModulus() -
                      n.Dot(c2) / c1.SquareModulus()) *
             input.scale;
    };
    // Never in the last refinements, which the positions keep.
    const bool g2 = input.curvature &&
                    bendRounds < input.bendRefinements && refinement < 8;
    // The error between the sites: section poles (homogeneous, a bound on the
    // section for equal weights), and the feet mapped back onto their supports.
    std::vector<double> refine;
    bool bent = false;
    double worstError = 0;
    for (int i = 0; i + 1 < count; ++i)
      for (double f : {0.25, 0.5, 0.75}) {
        const double u = sites[i] + (sites[i + 1] - sites[i]) * f;
        const auto got = evaluate(u), want = row(u);
        double error = 0;
        for (int j = 0; j < S; ++j) {
          const double wg = got[4 * j + 3], ww = want[4 * j + 3];
          if (!(wg > 0)) {
            error = HUGE_VAL;
            break;
          }
          const gp_Pnt pg(got[4 * j] / wg, got[4 * j + 1] / wg,
                          got[4 * j + 2] / wg);
          const gp_Pnt pw(want[4 * j] / ww, want[4 * j + 1] / ww,
                          want[4 * j + 2] / ww);
          error =
              std::max(error, pg.Distance(pw) + std::abs(wg - ww) * input.scale);
        }
        for (int side = 0; side < 2; ++side) {
          const gp_Pnt foot = input.contactPoint(side, got[4 * S + 2 * side],
                                     got[4 * S + 2 * side + 1]);
          const gp_Pnt footWant = input.contactPoint(side, want[4 * S + 2 * side],
                                         want[4 * S + 2 * side + 1]);
          error = std::max(error, foot.Distance(footWant));
        }
        worstError = std::max(worstError, error);
        bool split = error > tolerance;
        if (!split && g2)
          for (int side = 0; side < 2 && !split; ++side)
            if (!input.cliff[side] && alongBend(u, side) > input.bendTolerance)
              split = bent = true;
        if (split) {
          const auto right = knots.upper_bound(u);
          refine.push_back((std::prev(right)->first + right->first) / 2);
        }
      }
    // Refined for the bend alone, a bounded number of times: the positions
    // decide the fit.
    if (bent)
      ++bendRounds;
    if (refine.empty()) {
      // The surface and the curves, from the coefficients.
      std::vector<double> distinct;
      std::vector<int> multiplicity;
      for (const auto &[u, m] : knots) {
        distinct.push_back(u);
        multiplicity.push_back(m);
      }
      NCollection_Array1<double> uKnots(1, int(distinct.size()));
      NCollection_Array1<int> uMults(1, int(distinct.size()));
      for (size_t k = 0; k < distinct.size(); ++k) {
        uKnots(int(k) + 1) = distinct[k];
        uMults(int(k) + 1) = multiplicity[k];
      }
      NCollection_Array1<double> vKnots(1, int(reference.knots.size()));
      NCollection_Array1<int> vMults(1, int(reference.knots.size()));
      for (size_t k = 0; k < reference.knots.size(); ++k) {
        vKnots(int(k) + 1) = reference.knots[k];
        vMults(int(k) + 1) = reference.multiplicities[k];
      }
      NCollection_Array2<gp_Pnt> poles(1, count, 1, S);
      NCollection_Array2<double> weights(1, count, 1, S);
      bool rational = false;
      for (int i = 0; i < count; ++i)
        for (int j = 0; j < S; ++j) {
          const double w = coefficients[4 * j + 3][i];
          if (!(w > 0))
            fail("blend-fit", no_solution, "the blend's weights do not stay positive");
          rational |= std::abs(w - 1) > input.degenerate;
          poles(i + 1, j + 1) =
              gp_Pnt(coefficients[4 * j][i] / w, coefficients[4 * j + 1][i] / w,
                     coefficients[4 * j + 2][i] / w);
          weights(i + 1, j + 1) = w;
        }
      BlendFitResult out;
      out.blend =
          rational
              ? new Geom_BSplineSurface(poles, weights, uKnots, vKnots, uMults,
                                        vMults, kDegree, reference.degree)
              : new Geom_BSplineSurface(poles, uKnots, vKnots, uMults, vMults,
                                        kDegree, reference.degree);
      for (int side = 0; side < 2; ++side) {
        NCollection_Array1<gp_Pnt2d> contact(1, count);
        for (int i = 0; i < count; ++i)
          contact(i + 1) = gp_Pnt2d(coefficients[4 * S + 2 * side][i],
                                    coefficients[4 * S + 2 * side + 1][i]);
        out.contact[side] =
            new Geom2d_BSplineCurve(contact, uKnots, uMults, kDegree);
      }
      myResult=std::move(out); return true;
    }
    // Breaks the fit only halves its error towards (de Boor) are looked for once a smooth fit should
    // have converged: a foot's short excursion past a support's end between the march's stations, and
    // where the supplied section law changes form: C0 in the poles there.
    if (refinement == 2) {
        if (raise(input.discoveredBreaks()))
          return false;
    }
    if (refinement == 3 && !contactPilot && input.curvature) {
        const auto& regime = input.regime;
        constexpr int scan = 128;
        int before = regime(0);
        std::map<double, int> changes;
        for (int k = 1; k <= scan; ++k) {
            const int now = regime(double(k) / scan);
            if (now != before) {
                double low = double(k - 1) / scan, high = double(k) / scan;
                for (int i = 0; i < 60 && high - low > input.degenerate; ++i) {
                    const double middle = (low + high) / 2;
                    (regime(middle) == before ? low : high) = middle;
                }
                input.addBreak(changes, (low + high) / 2, kDegree);
            }
            before = now;
        }
        if (raise(changes))
          return false;
    }
    // Halving the spans takes a smooth blend's error down some 64 times each
    // (degree 5). An error still a thousand tolerances that does not halve in
    // three refinements is a jump in the data (a contact that folds); near the
    // tolerance a slow tail goes on.
    worst.push_back(worstError);
    if (worst.size() > 3 && worstError > 1e3 * tolerance &&
        worstError > worst[worst.size() - 4] / 2)
      { myStopped=true; return false; }
    for (double u : refine)
      knots.emplace(u, 1);
  return false;
}
BlendFitResult fitBlendSections(const BlendFitInput& input,const Message_ProgressRange& progress) {
 BlendSectionFitter fitter(input,input.knots);
 for(int refinement=0;refinement<16;++refinement) {
  checkCancelled(progress);
  if(fitter.Advance(refinement))return fitter.Result();
  if(fitter.IsStopped())break;
 }
 fail("blend-fit",no_solution,"the blend surface cannot meet the length tolerance");
}
}
