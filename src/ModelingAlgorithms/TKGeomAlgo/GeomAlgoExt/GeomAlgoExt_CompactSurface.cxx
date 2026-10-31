// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#include <GeomAlgoExt_CompactSurface.hxx>

#include <GeomAlgoExt_Execution.hxx>
#include <GeomAlgoExt_Interval.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>

namespace GeomAlgoExt {
namespace {
inline bool finite(const gp_Pnt& p) { return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()); }
// Limits affect only this optional optimization; the accepted original is retained.
constexpr int mostTrials = 512, mostRemovals = 16;
constexpr int64_t mostCoefficients = 262144, mostVisits = 8 * 1024 * 1024;

using GeomAlgoExt::interval::Range;
using GeomAlgoExt::interval::add;
using GeomAlgoExt::interval::subtract;
using GeomAlgoExt::interval::multiply;
using GeomAlgoExt::interval::valid;
using GeomAlgoExt::interval::roundUp;
using IntervalPole = std::array<Range, 3>;
struct IntervalNet {
    int degreeU, degreeV, nu, nv;
    std::vector<double> u, v;
    std::vector<IntervalPole> poles;
};
std::optional<IntervalNet> intervalNet(const occ::handle<Geom_BSplineSurface> &surface) {
    IntervalNet net{surface->UDegree(), surface->VDegree(), surface->NbUPoles(), surface->NbVPoles(), {}, {}, {}};
    for (int i=1;i<=surface->NbUKnots();++i)
        net.u.insert(net.u.end(),surface->UMultiplicity(i),surface->UKnot(i));
    for (int j=1;j<=surface->NbVKnots();++j)
        net.v.insert(net.v.end(),surface->VMultiplicity(j),surface->VKnot(j));
    if (net.u.size()!=static_cast<size_t>(net.nu+net.degreeU+1) ||
        net.v.size()!=static_cast<size_t>(net.nv+net.degreeV+1)) return {};
    for (int i=1;i<=net.nu;++i) for (int j=1;j<=net.nv;++j) {
        const gp_Pnt p=surface->Pole(i,j);
        if (!finite(p)) return {};
        net.poles.push_back({Range{p.X(),p.X()},Range{p.Y(),p.Y()},Range{p.Z(),p.Z()}});
    }
    return net;
}
// One Boehm insertion: Q_i=alpha_i P_i+(1-alpha_i)P_(i-1).
// All knot values start as exact stored-double intervals. Positive denominator and
// alpha intersection with [0,1] are certified before any coefficient is used.
bool insertInterval(const Message_ProgressRange &operation, IntervalNet &net, bool alongU, double knot, int64_t &visits) {
    auto &knots=alongU?net.u:net.v;
    const int degree=alongU?net.degreeU:net.degreeV, count=alongU?net.nu:net.nv;
    const int other=alongU?net.nv:net.nu;
    const int k=static_cast<int>(std::upper_bound(knots.begin(),knots.end(),knot)-knots.begin())-1;
    const int multiplicity=static_cast<int>(std::count(knots.begin(),knots.end(),knot));
    if (!(knot>knots[degree] && knot<knots[count]) || multiplicity>=degree ||
        k-degree<0 || k-multiplicity>=count) return false;
    const int64_t work=static_cast<int64_t>(other)*(degree-multiplicity);
    if (visits+work>mostVisits) return false;
    visits+=work;
    std::vector<Range> alpha;
    for (int i=k-degree+1;i<=k-multiplicity;++i) {
        const Range numerator=subtract({knot,knot},{knots[i],knots[i]});
        const Range denominator=subtract({knots[i+degree],knots[i+degree]},{knots[i],knots[i]});
        if (!valid(numerator)||!valid(denominator)||!(denominator.low>0)) return false;
        const Range ratio=GeomAlgoExt::interval::intersect(GeomAlgoExt::interval::dividePositive(numerator,denominator),{0,1});
        if (!valid(ratio)) return false;
        alpha.push_back(ratio);
    }
    std::vector<IntervalPole> next(static_cast<size_t>(count+1)*other);
    const auto old=[&](int i,int j)->const IntervalPole& {return net.poles[alongU?i*net.nv+j:j*net.nv+i];};
    const auto at=[&](int i,int j)->IntervalPole& {return next[alongU?i*net.nv+j:j*(net.nv+1)+i];};
    for (int j=0;j<other;++j) {
        if (j%16==0) checkCancelled(operation);
        for (int i=0;i<=k-degree;++i) at(i,j)=old(i,j);
        for (int i=k-multiplicity;i<count;++i) at(i+1,j)=old(i,j);
        for (int i=k-degree+1;i<=k-multiplicity;++i) {
            const Range a=alpha[i-(k-degree+1)], b=subtract({1,1},a);
            for (int axis=0;axis<3;++axis) {
                at(i,j)[axis]=add(multiply(a,old(i,j)[axis]),multiply(b,old(i-1,j)[axis]));
                if (!valid(at(i,j)[axis])) return false;
            }
        }
    }
    knots.insert(knots.begin()+k+1,knot);
    if (alongU) ++net.nu; else ++net.nv;
    net.poles=std::move(next);
    return true;
}
bool bezierIntervals(const Message_ProgressRange &operation, IntervalNet &net,
                     const occ::handle<Geom_BSplineSurface> &original, int64_t &visits) {
    for (int i=2;i<original->NbUKnots();++i)
        while (std::count(net.u.begin(),net.u.end(),original->UKnot(i))<net.degreeU)
            if (!insertInterval(operation,net,true,original->UKnot(i),visits)) return false;
    for (int j=2;j<original->NbVKnots();++j)
        while (std::count(net.v.begin(),net.v.end(),original->VKnot(j))<net.degreeV)
            if (!insertInterval(operation,net,false,original->VKnot(j),visits)) return false;
    return true;
}
// The work bezierIntervals counts for a net of these knots, and whether it would get through:
// the same insertions and the same checks of the knots, without the coefficients (which,
// finite, give finite intervals).
bool countIntervals(IntervalNet net, const occ::handle<Geom_BSplineSurface> &original, int64_t &visits) {
    for (const bool alongU : {true, false}) {
        auto &knots=alongU?net.u:net.v;
        const int degree=alongU?net.degreeU:net.degreeV;
        int &count=alongU?net.nu:net.nv;
        const int other=alongU?net.nv:net.nu;
        for (int i=2;i<(alongU?original->NbUKnots():original->NbVKnots());++i) {
            const double knot=alongU?original->UKnot(i):original->VKnot(i);
            while (std::count(knots.begin(),knots.end(),knot)<degree) {
                const int k=static_cast<int>(std::upper_bound(knots.begin(),knots.end(),knot)-knots.begin())-1;
                const int multiplicity=static_cast<int>(std::count(knots.begin(),knots.end(),knot));
                if (!(knot>knots[degree] && knot<knots[count]) || multiplicity>=degree ||
                    k-degree<0 || k-multiplicity>=count) return false;
                const int64_t work=static_cast<int64_t>(other)*(degree-multiplicity);
                if (visits+work>mostVisits) return false;
                visits+=work;
                for (int a=k-degree+1;a<=k-multiplicity;++a) {
                    const Range numerator=subtract({knot,knot},{knots[a],knots[a]});
                    const Range denominator=subtract({knots[a+degree],knots[a+degree]},{knots[a],knots[a]});
                    if (!valid(numerator)||!valid(denominator)||!(denominator.low>0)) return false;
                    if (!valid(GeomAlgoExt::interval::intersect(GeomAlgoExt::interval::dividePositive(numerator,denominator),{0,1}))) return false;
                }
                knots.insert(knots.begin()+k+1,knot);
                ++count;
            }
        }
    }
    return true;
}
// A distance that `trial` certainly departs from `original` by somewhere: their points apart
// along the line of the knot removed, less what rounding of the two evaluations can make of
// it. The common-basis bound is never below the surfaces' distance at any parameters.
double departure(const occ::handle<Geom_BSplineSurface> &original, const occ::handle<Geom_BSplineSurface> &trial,
                 bool alongU, double knot) {
    double largest = 0;
    for (const auto *surface : {&original, &trial})
        for (int i = 1; i <= (*surface)->NbUPoles(); ++i)
            for (int j = 1; j <= (*surface)->NbVPoles(); ++j) {
                const gp_Pnt p = (*surface)->Pole(i, j);
                largest = std::max({largest, std::abs(p.X()), std::abs(p.Y()), std::abs(p.Z())});
            }
    double u0, u1, v0, v1;
    original->Bounds(u0, u1, v0, v1);
    double apart = 0;
    for (int k = 0; k <= 32; ++k) {
        const double u = alongU ? knot : u0 + (u1 - u0) * k / 32, v = alongU ? v0 + (v1 - v0) * k / 32 : knot;
        apart = std::max(apart, original->Value(u, v).Distance(trial->Value(u, v)));
    }
    return apart - 1e-10 * (1 + largest);
}

double intervalDisplacement(const Message_ProgressRange &operation, const IntervalNet &reference, const IntervalNet &candidate) {
    // Flat knot equality checks every unique knot and its multiplicity, not only dimensions.
    if (reference.degreeU!=candidate.degreeU || reference.degreeV!=candidate.degreeV ||
        reference.nu!=candidate.nu || reference.nv!=candidate.nv ||
        reference.u!=candidate.u || reference.v!=candidate.v) return std::numeric_limits<double>::infinity();
    double worst=0;
    for (size_t i=0;i<reference.poles.size();++i) {
        if (i%256==0) checkCancelled(operation);
        double squared=0;
        for (int axis=0;axis<3;++axis) {
            const Range difference=subtract(reference.poles[i][axis],candidate.poles[i][axis]);
            if (!valid(difference)) return std::numeric_limits<double>::infinity();
            const double magnitude=std::max(std::abs(difference.low),std::abs(difference.high));
            squared=roundUp(squared+roundUp(magnitude*magnitude));
        }
        const double distance=roundUp(std::sqrt(squared));
        if (!std::isfinite(distance)) return std::numeric_limits<double>::infinity();
        worst=std::max(worst,distance);
    }
    return worst;
}

}

std::optional<CompactPatchSurface> compactSurface(
    const Message_ProgressRange &operation, const occ::handle<Geom_BSplineSurface> &original,
    double maximumDisplacement) {
    if (original.IsNull() || !(maximumDisplacement > 0) || !std::isfinite(maximumDisplacement) ||
        original->IsURational() || original->IsVRational() ||
        original->IsUPeriodic() || original->IsVPeriodic() ||
        original->UMultiplicity(1)!=original->UDegree()+1 ||
        original->UMultiplicity(original->NbUKnots())!=original->UDegree()+1 ||
        original->VMultiplicity(1)!=original->VDegree()+1 ||
        original->VMultiplicity(original->NbVKnots())!=original->VDegree()+1) return {};
    const int64_t rows = 1 + static_cast<int64_t>(original->NbUKnots() - 1) * original->UDegree();
    const int64_t columns = 1 + static_cast<int64_t>(original->NbVKnots() - 1) * original->VDegree();
    if (rows < 2 || columns < 2 || rows > mostCoefficients / columns) return {};
    for (int i = 1; i <= original->NbUKnots(); ++i)
        if (!std::isfinite(original->UKnot(i)) || (i>1 && !(original->UKnot(i)>original->UKnot(i-1))) ||
            original->UMultiplicity(i)<1 || (i>1 && i<original->NbUKnots() && original->UMultiplicity(i)>original->UDegree())) return {};
    for (int j = 1; j <= original->NbVKnots(); ++j)
        if (!std::isfinite(original->VKnot(j)) || (j>1 && !(original->VKnot(j)>original->VKnot(j-1))) ||
            original->VMultiplicity(j)<1 || (j>1 && j<original->NbVKnots() && original->VMultiplicity(j)>original->VDegree())) return {};
    double diameter = 0;
    const gp_Pnt anchor = original->Pole(1, 1);
    for (int i = 1; i <= original->NbUPoles(); ++i) {
        if (i % 16 == 0) checkCancelled(operation);
        for (int j = 1; j <= original->NbVPoles(); ++j) {
            const gp_Pnt p = original->Pole(i, j);
            if (!finite(p)) return {};
            diameter = std::max(diameter, p.Distance(anchor));
        }
    }
    // OCCT's removal estimate only proposes a surface. The independent common-basis
    // bound below, not this deliberately loose estimate, controls physical displacement.
    const double proposalTolerance = 1024 * std::max(diameter, maximumDisplacement);
    if (!std::isfinite(proposalTolerance)) return {};
    try {
        int64_t visits = 0;
        auto reference = intervalNet(original);
        if (!reference || !bezierIntervals(operation,*reference,original,visits)) return {};
        const int64_t coefficients = static_cast<int64_t>(reference->nu) * reference->nv;
        if (coefficients > mostCoefficients || coefficients < 1) return {};
        auto current = occ::down_cast<Geom_BSplineSurface>(original->Copy());
        int trials = 0, removals = 0;
        double displacement = 0;
        while (removals < mostRemovals && trials < mostTrials && visits + coefficients <= mostVisits) {
            checkCancelled(operation);
            occ::handle<Geom_BSplineSurface> best;
            double least = std::numeric_limits<double>::infinity();
            for (int axis = 0; axis < 2; ++axis) {
                const int knots = axis == 0 ? current->NbUKnots() : current->NbVKnots();
                for (int knot = 2; knot < knots && trials < mostTrials && visits + coefficients <= mostVisits; ++knot) {
                    checkCancelled(operation);
                    ++trials;
                    auto trial = occ::down_cast<Geom_BSplineSurface>(current->Copy());
                    try {
                        if (!(axis == 0 ? trial->RemoveUKnot(knot, 0, proposalTolerance)
                                        : trial->RemoveVKnot(knot, 0, proposalTolerance))) continue;
                        auto expanded = intervalNet(trial);
                        if (!expanded) continue;
                        // Too far from the original where it is read: refused as the common basis
                        // would refuse it, with the work it would have counted.
                        if (departure(original, trial, axis == 0,
                                      axis == 0 ? current->UKnot(knot) : current->VKnot(knot)) > maximumDisplacement) {
                            if (countIntervals({expanded->degreeU, expanded->degreeV, expanded->nu, expanded->nv,
                                                expanded->u, expanded->v, {}}, original, visits) &&
                                visits + coefficients <= mostVisits)
                                visits += coefficients;
                            continue;
                        }
                        if (!bezierIntervals(operation,*expanded,original,visits)) continue;
                        if (visits+coefficients>mostVisits) continue;
                        visits += coefficients;
                        const double error = intervalDisplacement(operation,*reference,*expanded);
                        if (error <= maximumDisplacement && error < least) { best = trial; least = error; }
                    } catch (const GeomAlgoExt_Failure &) { throw; } catch (const Standard_Failure &) { /* Keep the original on an unavailable proposal. */ }
                }
            }
            if (best.IsNull()) break;
            current = best;
            displacement = least;
            ++removals;
        }
        if (removals) return CompactPatchSurface{current, displacement, removals};
    } catch (const GeomAlgoExt_Failure &) { throw; } catch (const Standard_Failure &) { /* Optional optimization cannot invalidate its input. */ }
    return {};
}
}
