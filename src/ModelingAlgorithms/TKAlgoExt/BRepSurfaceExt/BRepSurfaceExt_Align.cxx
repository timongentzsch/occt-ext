// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

namespace loft {
EvolutionData alignSurface(OcctKernel &k, uint32_t faceId, uint32_t edgeId, uint32_t referenceId,
                           uint32_t supportId, int order, double tolerance, int extraRows, int pu, int pv, int su, int sv, bool normalFlow,
                           double partialStart, double partialEnd, bool skewed, double skewStart, double skewEnd) {
    const char *op = "alignSurface";
    if (!std::isfinite(tolerance) || tolerance <= 0 || order < 0 || order > 2)
        fail(op, invalid_input, "continuity must be G0/G1/G2 and tolerance positive and finite");
    if (!(partialStart >= 0 && partialStart < partialEnd && partialEnd <= 1))
        fail(op, invalid_input, "the partial interval must rise within 0 to 1");
    if (skewed && order < 1)
        fail(op, invalid_input, "a skewed boundary needs at least G1");
    if (skewed && !(std::abs(skewStart) < M_PI / 2 && std::abs(skewEnd) < M_PI / 2))
        fail(op, invalid_input, "skew angles must lie within a right angle of straight across");
    if (extraRows < 0 || extraRows > 9 || pu < 0 || pu > 5 || pv < 0 || pv > 5 || su < 0 || su > 9 || sv < 0 || sv > 9)
        fail(op, invalid_input, "extra rows 0..9, optional degrees 1..5 and spans 1..9 are required");
    const auto face = surfacing::faceOf(k, op, faceId);
    const auto own = surfacing::boundaryOf(k, op, edgeId, faceId);
    const auto ref = surfacing::boundaryOf(k, op, referenceId, supportId);
    if (face.IsSame(ref.face))
        fail(op, invalid_input, "the reference must be a different face");
    auto surface = surfacing::editable(op, face);
    if ((pu && pu < surface->UDegree()) || (pv && pv < surface->VDegree()))
        fail(op, unsupported_geometry, "Align can raise degree exactly; lower it with Rebuild Surface first");
    surface->IncreaseDegree(pu ? pu : surface->UDegree(), pv ? pv : surface->VDegree());
    // Knot insertion is exact. Split the widest current interval until the requested
    // span count is reached; never silently discard existing breakpoints.
    const auto refine = [&](bool u, int spans) {
        const auto count = [&] { return (u ? surface->NbUKnots() : surface->NbVKnots()) - 1; };
        if (spans && spans < count())
            fail(op, unsupported_geometry, "Align can add spans exactly; reduce them with Rebuild Surface first");
        while (count() < spans) {
            int at = 1;
            const auto knot = [&](int i) { return u ? surface->UKnot(i) : surface->VKnot(i); };
            for (int i = 2; i <= count(); ++i)
                if (knot(i + 1) - knot(i) > knot(at + 1) - knot(at)) at = i;
            const double middle = (knot(at) + knot(at + 1)) / 2;
            if (u) surface->InsertUKnot(middle, 1, 1e-12); else surface->InsertVKnot(middle, 1, 1e-12);
        }
    };
    refine(true, su); refine(false, sv);
    auto net = surfacing::netOf(surface);
    if (net.poles.size() > 1024)
        fail(op, unsupported_geometry, "local fitting is limited to 1024 poles");
    const int side = surfacing::boundarySide(op, face, own.edge), rows = order + 1 + extraRows;
    const int nu = net.countU(), nv = net.countV();
    if (extraRows == 0 && (side < 2 ? nu : nv) <= rows)
        fail(op, unsupported_geometry,
             "not enough rows to match while keeping the opposite boundary fixed; raise degree first");
    std::vector<char> fixed(net.poles.size(), 1);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            const int from = side == 0 ? i : side == 1 ? nu - 1 - i : side == 2 ? j : nv - 1 - j;
            fixed[i * nv + j] = from >= rows;
        }
    auto target = surfacing::constraint(op, face, own.edge, ref, order);
    // Only this stretch of the reference, its shares counted from the start of the moving edge.
    const double whole = target.reference.last - target.reference.first;
    target.reference.last = target.reference.first + whole * partialEnd;
    target.reference.first = target.reference.first + whole * partialStart;
    target.skewed = skewed;
    target.skewStart = skewed ? std::tan(skewStart) : 0;
    target.skewEnd = skewed ? std::tan(skewEnd) : 0;
    const std::vector<surfacing::Constraint> targets{target};
    // Exact degree/knot refinement should not turn an already matched surface into a
    // different one merely to reduce the fitting energy.
    try {
        const auto unchanged = surfacing::makeFace(op, surface, tolerance, face.Orientation());
        surfacing::measure(op, unchanged, targets, tolerance, normalFlow);
        return stored(k, unchanged);
    } catch (const Failure &failure) {
        if (failure.reason != no_solution) throw;
    }
    for (int pass = 0; pass < (order == 2 || skewed ? 8 : 1); ++pass) {
        fitNet(op, net, fixed, surfacing::conditions(op, targets, surface, tolerance, normalFlow), {1e-11, 1e-13}, {});
        surface = surfaceOf(net);
    }
    const auto made = surfacing::makeFace(op, surface, tolerance, face.Orientation());
    surfacing::measure(op, made, targets, tolerance, normalFlow);
    return stored(k, made);
}

namespace {
/// One edge of a spline face, as what the Collinear solver needs of it: which side, and the poles
/// of the boundary row and the next one after they are brought to a common basis.
struct SeamSide {
    occ::handle<Geom_BSplineSurface> surface;
    int side;
    bool alongV() const { return side < 2; }
    int rows() const { return alongV() ? surface->NbVPoles() : surface->NbUPoles(); }
    int depth(int k) const { return side == 0 || side == 2 ? k : (alongV() ? surface->NbUPoles() : surface->NbVPoles()) + 1 - k; }
    gp_Pnt pole(int row, int k) const { return alongV() ? surface->Pole(depth(row), k) : surface->Pole(k, depth(row)); }
    void setPole(int row, int k, const gp_Pnt &p) {
        if (alongV())
            surface->SetPole(depth(row), k, p);
        else
            surface->SetPole(k, depth(row), p);
    }
    std::pair<double, double> range() const {
        double a, b, c, d;
        surface->Bounds(a, b, c, d);
        return alongV() ? std::pair{c, d} : std::pair{a, b};
    }
    int degree() const { return alongV() ? surface->VDegree() : surface->UDegree(); }
    int knots() const { return alongV() ? surface->NbVKnots() : surface->NbUKnots(); }
    double knot(int k) const { return alongV() ? surface->VKnot(k) : surface->UKnot(k); }
    int multiplicity(int k) const { return alongV() ? surface->VMultiplicity(k) : surface->UMultiplicity(k); }
    double share(int k) const {
        const auto [a, b] = range();
        return (knot(k) - a) / (b - a);
    }
    void insert(double s, int multiplicity) {
        const auto [a, b] = range();
        if (alongV())
            surface->InsertVKnot(a + s * (b - a), multiplicity, 1e-12, false);
        else
            surface->InsertUKnot(a + s * (b - a), multiplicity, 1e-12, false);
    }
    void raise(int to) {
        if (alongV())
            surface->IncreaseDegree(surface->UDegree(), to);
        else
            surface->IncreaseDegree(to, surface->VDegree());
    }
};
} // namespace

// The Collinear solver of Edge Align: along the edge shared with a reference, each pole of the row
// next to the boundary is put on the line through the boundary pole and the reference's next pole,
// so the hulls of the two surfaces make one grid. Asymmetric, only the moving surface changes:
// its boundary on the reference's and its next row in line with the reference's. Symmetric, both
// change: the boundary halfway, and each row's lines along their shared mean direction. Each next
// pole keeps its distance from the boundary, or with `matchRatio` all have the one ratio to the
// reference's, which makes the join exactly tangent. A skew leans the lines along the
// boundary, as in the exact solve. What it gives is one face, or with `symmetric` a compound of
// the moving face and the reference.
EvolutionData alignCollinear(OcctKernel &k, uint32_t faceId, uint32_t edgeId, uint32_t referenceId, uint32_t supportId,
                             bool symmetric, bool matchRatio, double partialStart, double partialEnd, double skewStart,
                             double skewEnd, double tolerance) {
    const char *op = "alignCollinear";
    if (!std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "tolerance must be positive and finite");
    if (!(partialStart >= 0 && partialStart < partialEnd && partialEnd <= 1))
        fail(op, invalid_input, "the partial interval must rise within 0 to 1");
    if (!(std::abs(skewStart) < M_PI / 2 && std::abs(skewEnd) < M_PI / 2))
        fail(op, invalid_input, "skew angles must lie within a right angle of straight across");
    const bool partial = partialStart > 0 || partialEnd < 1;
    if (symmetric && partial)
        fail(op, invalid_input, "only the moving surface can follow a part of the reference; choose Asymmetric");
    const auto faceShape = surfacing::faceOf(k, op, faceId);
    const auto own = surfacing::boundaryOf(k, op, edgeId, faceId);
    const auto ref = surfacing::boundaryOf(k, op, referenceId, supportId);
    if (faceShape.IsSame(ref.face))
        fail(op, invalid_input, "the reference must be a different face");
    SeamSide moving{surfacing::editable(op, faceShape), surfacing::boundarySide(op, faceShape, own.edge)};
    SeamSide fixed{surfacing::editable(op, ref.face), surfacing::boundarySide(op, ref.face, ref.edge)};
    const auto endPoint = [](const SeamSide &s, bool last) {
        const auto [a, b] = s.range();
        double u0, u1, v0, v1;
        s.surface->Bounds(u0, u1, v0, v1);
        const double across = s.side == 0 ? u0 : s.side == 1 ? u1 : s.side == 2 ? v0 : v1;
        const double along = last ? b : a;
        return s.alongV() ? s.surface->Value(across, along) : s.surface->Value(along, across);
    };
    // Whether the reference runs the other way along the edge.
    const gp_Pnt a0 = endPoint(moving, false), a1 = endPoint(moving, true), b0 = endPoint(fixed, false), b1 = endPoint(fixed, true);
    const bool opposed = a0.Distance(b1) + a1.Distance(b0) < a0.Distance(b0) + a1.Distance(b1);
    // A part of the reference is the whole of a copy: the moving row is laid on it.
    if (partial) {
        const auto [a, b] = fixed.range();
        const double from = opposed ? 1 - partialEnd : partialStart, to = opposed ? 1 - partialStart : partialEnd;
        double u0, u1, v0, v1;
        fixed.surface->Bounds(u0, u1, v0, v1);
        if (fixed.alongV())
            fixed.surface->Segment(u0, u1, a + from * (b - a), a + to * (b - a));
        else
            fixed.surface->Segment(a + from * (b - a), a + to * (b - a), v0, v1);
    }
    // One degree and one set of knots along the edge, by exact refinement of both.
    const int degree = std::max(moving.degree(), fixed.degree());
    moving.raise(degree);
    fixed.raise(degree);
    std::vector<std::pair<double, int>> movingKnots, fixedKnots;
    for (int i = 2; i < moving.knots(); ++i)
        movingKnots.emplace_back(moving.share(i), moving.multiplicity(i));
    for (int i = 2; i < fixed.knots(); ++i)
        fixedKnots.emplace_back(fixed.share(i), fixed.multiplicity(i));
    for (const auto &[s, m] : movingKnots)
        fixed.insert(opposed ? 1 - s : s, m);
    for (const auto &[s, m] : fixedKnots)
        moving.insert(opposed ? 1 - s : s, m);
    const int count = moving.rows();
    if (count != fixed.rows())
        fail(op, no_solution, "the two surfaces could not be brought to one basis along the edge");
    const auto partner = [&](int i) { return opposed ? count + 1 - i : i; };
    double ratio = 0;
    if (matchRatio) {
        for (int i = 1; i <= count; ++i)
            ratio += moving.pole(2, i).Distance(moving.pole(1, i)) / fixed.pole(2, partner(i)).Distance(fixed.pole(1, partner(i))) / count;
    }
    const double leanStart = std::tan(skewStart), leanEnd = std::tan(skewEnd);
    for (int i = 1; i <= count; ++i) {
        const gp_Pnt m0 = moving.pole(1, i), m1 = moving.pole(2, i), f0 = fixed.pole(1, partner(i)), f1 = fixed.pole(2, partner(i));
        const double movingLength = m0.Distance(m1), fixedLength = f0.Distance(f1);
        if (movingLength <= gp::Resolution() || fixedLength <= gp::Resolution())
            fail(op, no_solution, "a row next to the edge lies on the edge");
        const gp_Vec fromFixed(f1, f0);
        const gp_Vec fromMoving(m0, m1);
        gp_Vec way = symmetric ? fromFixed.Normalized() + fromMoving.Normalized() : fromFixed;
        if (way.Magnitude() <= 1e-9 * fixedLength)
            fail(op, no_solution, "the two surfaces fold back on each other along the edge");
        way.Normalize();
        const gp_Pnt seam = symmetric ? gp_Pnt((m0.XYZ() + f0.XYZ()) / 2) : f0;
        const double length = matchRatio ? ratio * fixedLength : movingLength;
        // The lean: along the boundary, towards the other end when positive.
        const double f = (i - 1.0) / std::max(count - 1, 1);
        const double lean = leanStart * (1 - f) - leanEnd * f;
        gp_Vec running(moving.pole(1, std::max(1, i - 1)), moving.pole(1, std::min(count, i + 1)));
        gp_Vec leaning = way;
        if (lean != 0 && running.Magnitude() > gp::Resolution()) {
            running.Normalize();
            running -= way * running.Dot(way);
            leaning = (way + running.Normalized() * lean).Normalized();
        }
        moving.setPole(1, i, seam);
        moving.setPole(2, i, seam.Translated(leaning * length));
        if (symmetric) {
            fixed.setPole(1, partner(i), seam);
            fixed.setPole(2, partner(i), seam.Translated(way * -fixedLength));
        }
    }
    const auto made = surfacing::makeFace(op, moving.surface, tolerance, faceShape.Orientation());
    if (!symmetric)
        return stored(k, made);
    TopoDS_Compound both;
    BRep_Builder builder;
    builder.MakeCompound(both);
    builder.Add(both, made);
    builder.Add(both, surfacing::makeFace(op, fixed.surface, tolerance, ref.face.Orientation()));
    return stored(k, both);
}
} // namespace loft
