// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Fillet v2: Where curves and surfaces meet, the edges made there, and the body's edges and faces cut back to them.

#include <BRepBlendExt_Internal.hxx>

#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRep_TVertex.hxx>
#include <GProp_GProps.hxx>

#include <functional>

namespace loft::blend2::detail {

TopoDS_Edge trimmed(const TopoDS_Edge& edge, const TopoDS_Vertex& from, const TopoDS_Vertex& to, double at) {
    const TopoDS_Edge forward = TopoDS::Edge(edge.Oriented(TopAbs_FORWARD));
    TopoDS_Vertex first, last;
    TopExp::Vertices(forward, first, last);
    double f, l;
    BRep_Tool::Range(forward, f, l);
    const bool atFirst = first.IsSame(from);
    BRep_Builder builder;
    TopoDS_Edge out = TopoDS::Edge(forward.EmptyCopied());
    builder.Add(out, (atFirst ? to : first).Oriented(TopAbs_FORWARD));
    builder.Add(out, (atFirst ? last : to).Oriented(TopAbs_REVERSED));
    builder.Range(out, atFirst ? at : f, atFirst ? l : at);
    builder.UpdateVertex(to, at, out, BRep_Tool::Tolerance(to));
    return out;
}

std::vector<TopoDS_Face> rebuiltFaces(const TopoDS_Face& face, const EdgeMap& replace, const std::vector<TopoDS_Edge>& bridges) {
    BRep_Builder builder;
    TopoDS_Face out = TopoDS::Face(face.EmptyCopied());
    // Wires and edges as the face stores them (their own orientations, where they are placed),
    // so that each comes back stored as it was, as OpenCascade's faces keep theirs.
    const auto stored = [&](const TopoDS_Shape& wire) { return wire.Oriented(TopAbs::Compose(wire.Orientation(), out.Orientation())); };
    // Each wire that changes, its edges as they run in the face; the bridges join each to the wire
    // whose gap they close, and two wires into one where they close a gap between them (a blend that
    // takes the strip between a hole and the face's boundary whole).
    std::vector<std::vector<TopoDS_Edge>> wires;
    std::vector<int> parent;
    const std::function<int(int)> root = [&](int i) { return parent[i] == i ? i : parent[i] = root(parent[i]); };
    NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> wireOf;
    // How many more edges arrive at each vertex than leave it: zero round closed wires.
    NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> balance;
    const auto add = [&](int w, const TopoDS_Shape& edge) {
        const TopoDS_Edge e = TopoDS::Edge(edge);
        wires[w].push_back(e);
        for (const auto& [vertex, step] : {std::pair{TopExp::LastVertex(e, true), 1}, std::pair{TopExp::FirstVertex(e, true), -1}}) {
            if (!balance.IsBound(vertex)) balance.Bind(vertex, 0);
            balance.ChangeFind(vertex) += step;
            if (!wireOf.IsBound(vertex)) wireOf.Bind(vertex, w);
            else if (root(wireOf.Find(vertex)) != root(w)) parent[root(wireOf.Find(vertex))] = root(w);
        }
    };
    for (TopoDS_Iterator w(face, false, true); w.More(); w.Next()) {
        bool touched = false;
        for (TopoDS_Iterator e(w.Value(), false, true); e.More(); e.Next()) touched |= replace.IsBound(e.Value());
        if (!touched) {
            builder.Add(out, stored(w.Value()));
            continue;
        }
        const int index = int(wires.size());
        wires.emplace_back();
        parent.push_back(index);
        const TopAbs_Orientation orientation = w.Value().Orientation();
        for (TopoDS_Iterator e(w.Value(), false, true); e.More(); e.Next()) {
            // An edge replaced by nothing is gone: the bridges close where it was.
            if (replace.IsBound(e.Value()) && replace.Find(e.Value()).IsNull()) continue;
            // A replacement that runs against the edge it replaces is bound reversed.
            const TopoDS_Shape edge = !replace.IsBound(e.Value())
                                          ? e.Value()
                                          : replace.Find(e.Value()).Oriented(TopAbs::Compose(e.Value().Orientation(), replace.Find(e.Value()).Orientation()));
            // In parts where it is bound to a compound: a seam left in two.
            if (edge.ShapeType() == TopAbs_COMPOUND) {
                for (TopoDS_Iterator part(edge); part.More(); part.Next())
                    add(index, part.Value().Oriented(TopAbs::Compose(e.Value().Orientation(), orientation)));
                continue;
            }
            add(index, edge.Oriented(TopAbs::Compose(edge.Orientation(), orientation)));
        }
    }
    // The bridges close the gaps the replaced edges leave: each runs on from a vertex that more
    // edges arrive at than leave, so that the wire runs one way round, whichever order its edges
    // came in (at a vertex the wire passes twice, a closed edge's seam, as at any).
    const auto surplus = [&](const TopoDS_Vertex& v) { return balance.IsBound(v) && balance.Find(v) > 0; };
    std::vector<bool> used(bridges.size(), false);
    for (bool progress = true; progress;) {
        progress = false;
        for (size_t i = 0; i < bridges.size(); ++i) {
            if (used[i]) continue;
            TopoDS_Vertex a, z;
            TopExp::Vertices(bridges[i], a, z);
            if (!surplus(a) && !surplus(z)) continue;
            const TopoDS_Vertex& from = surplus(a) ? a : z;
            if (!wireOf.IsBound(from)) continue;
            used[i] = progress = true;
            add(root(wireOf.Find(from)), bridges[i].Oriented(surplus(a) ? TopAbs_FORWARD : TopAbs_REVERSED));
        }
    }
    // Bridges left over that close among themselves, round places where nothing of the wire they
    // replace is left (the islands a neck that the blends close leaves of a cap), are wires of their
    // own; a bridge that closes on itself is kept with the first wire.
    std::vector<bool> island(wires.size(), false);
    for (size_t i = 0; i < bridges.size(); ++i) {
        if (used[i] || TopExp::FirstVertex(bridges[i]).IsSame(TopExp::LastVertex(bridges[i]))) continue;
        std::vector<TopoDS_Edge> loop = {bridges[i]};
        std::vector<size_t> taken = {i};
        TopoDS_Vertex start = TopExp::FirstVertex(bridges[i]), at = TopExp::LastVertex(bridges[i]);
        while (!at.IsSame(start)) {
            size_t next = bridges.size();
            for (size_t j = 0; j < bridges.size() && next == bridges.size(); ++j)
                if (!used[j] && std::find(taken.begin(), taken.end(), j) == taken.end() &&
                    (TopExp::FirstVertex(bridges[j]).IsSame(at) || TopExp::LastVertex(bridges[j]).IsSame(at)))
                    next = j;
            if (next == bridges.size()) break;
            const bool forward = TopExp::FirstVertex(bridges[next]).IsSame(at);
            loop.push_back(TopoDS::Edge(bridges[next].Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED)));
            taken.push_back(next);
            at = forward ? TopExp::LastVertex(bridges[next]) : TopExp::FirstVertex(bridges[next]);
        }
        if (!at.IsSame(start)) continue;
        for (size_t j : taken) used[j] = true;
        // Each edge's curve on the face runs on from the one before, in its period: pieces of a chain
        // whose vertex lay on a seam that has gone start a period apart.
        for (size_t j = 1; j < loop.size(); ++j) {
            double f, l, g, h;
            const auto last = BRep_Tool::CurveOnSurface(loop[j - 1], out, f, l), next = BRep_Tool::CurveOnSurface(loop[j], out, g, h);
            if (last.IsNull() || next.IsNull()) continue;
            const gp_Pnt2d end = last->Value(loop[j - 1].Orientation() == TopAbs_FORWARD ? l : f);
            const auto laid = aligned(next, loop[j].Orientation() == TopAbs_FORWARD ? g : h, end, BRep_Tool::Surface(out));
            if (laid != next) builder.UpdateEdge(loop[j], laid, out, BRep_Tool::Tolerance(loop[j]));
        }
        wires.push_back(loop);
        parent.push_back(int(parent.size()));
        island.push_back(true);
    }
    std::vector<TopoDS_Shape> islands;
    bool first = true;
    for (size_t w = 0; w < wires.size(); ++w) {
        if (root(int(w)) != int(w)) continue;
        TopoDS_Wire wire;
        builder.MakeWire(wire);
        int count = 0;
        for (size_t m = 0; m < wires.size(); ++m)
            if (root(int(m)) == int(w))
                for (const auto& edge : wires[m]) {
                    builder.Add(wire, edge);
                    ++count;
                }
        // A wire whose edges are all gone, its place taken by islands of bridges.
        if (count == 0) continue;
        if (first)
            for (size_t i = 0; i < bridges.size(); ++i)
                if (!used[i] && TopExp::FirstVertex(bridges[i]).IsSame(TopExp::LastVertex(bridges[i]))) builder.Add(wire, bridges[i]);
        first = false;
        // The parts that share no vertex (the blend took the strip between them, the contact out
        // over the face's edge and back) are islands of their own.
        std::vector<TopoDS_Edge> left;
        for (TopoDS_Iterator e(wire, false, true); e.More(); e.Next()) left.push_back(TopoDS::Edge(e.Value()));
        std::vector<std::vector<TopoDS_Edge>> parts;
        while (!left.empty()) {
            std::vector<TopoDS_Edge> part = {left.back()};
            left.pop_back();
            for (bool grew = true; grew;) {
                grew = false;
                for (size_t i = 0; i < left.size();) {
                    const bool touches = std::any_of(part.begin(), part.end(), [&](const TopoDS_Edge& e) {
                        return TopExp::FirstVertex(e).IsSame(TopExp::FirstVertex(left[i])) || TopExp::FirstVertex(e).IsSame(TopExp::LastVertex(left[i])) ||
                               TopExp::LastVertex(e).IsSame(TopExp::FirstVertex(left[i])) || TopExp::LastVertex(e).IsSame(TopExp::LastVertex(left[i]));
                    });
                    if (!touches) {
                        ++i;
                        continue;
                    }
                    part.push_back(left[i]);
                    left.erase(left.begin() + i);
                    grew = true;
                }
            }
            parts.push_back(part);
        }
        for (const auto& part : parts) {
            TopoDS_Wire one = wire;
            if (parts.size() > 1) {
                builder.MakeWire(one);
                for (const auto& edge : part) builder.Add(one, edge);
            }
            one.Closed(true);
            builder.Add(out, stored(one));
            if (island[w] || parts.size() > 1) islands.push_back(one);
        }
    }
    // Each wire closes round its edges and each edge lies on the face: a face that does not is
    // refused here, before OpenCascade's own checks read it (which, with their guards compiled out,
    // would follow its missing curves through null).
    for (TopoDS_Iterator w(out); w.More(); w.Next()) {
        NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> through;
        for (TopExp_Explorer e(w.Value(), TopAbs_EDGE); e.More(); e.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(e.Current());
            double f, l;
            if (!BRep_Tool::Degenerated(edge) && BRep_Tool::CurveOnSurface(edge, out, f, l).IsNull())
                throw Failure("filletV2", no_solution, "an edge of a face the blend changes does not lie on it", {});
            for (const auto& [vertex, step] : {std::pair{TopExp::LastVertex(edge, true), 1}, std::pair{TopExp::FirstVertex(edge, true), -1}}) {
                if (vertex.IsNull()) throw Failure("filletV2", no_solution, "an edge of a face the blend changes has no vertex at an end", {});
                if (!through.IsBound(vertex)) through.Bind(vertex, 0);
                through.ChangeFind(vertex) += step;
            }
        }
        for (NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher>::Iterator it(through); it.More(); it.Next())
            if (it.Value() != 0) throw Failure("filletV2", no_solution, "a face the blend changes does not close round its edges", {});
    }
    if (islands.empty()) return {out};
    // Islands: a face for each wire that bounds from outside (its own area positive), with the
    // wires inside it. An island's wire, which its bridges alone make, runs the way that bounds it
    // from outside.
    const auto alone = [&](const TopoDS_Shape& wire) {
        TopoDS_Face one = TopoDS::Face(out.EmptyCopied());
        builder.Add(one, wire);
        return one;
    };
    std::vector<TopoDS_Face> faces;
    std::vector<TopoDS_Shape> inner;
    for (TopoDS_Iterator w(out); w.More(); w.Next()) {
        GProp_GProps area;
        BRepGProp::SurfaceProperties(alone(w.Value()), area);
        const bool own = std::any_of(islands.begin(), islands.end(), [&](const TopoDS_Shape& i) { return i.IsSame(w.Value()); });
        if (area.Mass() > 0) faces.push_back(alone(w.Value()));
        else if (own) faces.push_back(alone(w.Value().Reversed()));
        else inner.push_back(w.Value());
    }
    if (faces.empty()) throw Failure("filletV2", no_solution, "a face the blend changes has no outer boundary", {});
    for (const auto& wire : inner) {
        // The face whose outer boundary holds a place of the wire.
        const TopoDS_Edge edge = TopoDS::Edge(TopExp_Explorer(wire, TopAbs_EDGE).Current());
        double f, l;
        const auto pcurve = BRep_Tool::CurveOnSurface(edge, out, f, l);
        const auto holder = std::find_if(faces.begin(), faces.end(), [&](const TopoDS_Face& candidate) {
            return BRepClass_FaceClassifier(candidate, pcurve->Value((f + l) / 2), tol::length).State() != TopAbs_OUT;
        });
        if (holder == faces.end()) throw Failure("filletV2", no_solution, "a hole of a face the blend changes lies in none of its islands", {});
        builder.Add(*holder, wire);
    }
    return faces;
}

bool poleInStrip(const Piece& piece, int side, const TopoDS_Edge& seam) {
    const TopoDS_Face& face = piece.rail.s[side].face;
    gp_Pnt pole;
    bool found = false;
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More() && !found; e.Next())
        if (BRep_Tool::Degenerated(TopoDS::Edge(e.Current())))
            for (TopExp_Explorer v(seam, TopAbs_VERTEX); v.More() && !found; v.Next())
                if (TopExp::FirstVertex(TopoDS::Edge(e.Current())).IsSame(v.Current())) {
                    pole = BRep_Tool::Pnt(TopoDS::Vertex(v.Current()));
                    found = true;
                }
    if (!found) return false;
    constexpr int samples = 256;
    const double u0 = piece.uAt(0), u1 = piece.uAt(1);
    double least = HUGE_VAL, nearest = u0;
    gp_Pnt foot;
    for (int i = 0; i <= samples; ++i) {
        const double u = u0 + (u1 - u0) * i / samples;
        const gp_Pnt2d uv = piece.fitted.contact[side]->Value(u);
        const gp_Pnt at = supportPoint(piece.rail.s[side].surface, uv.X(), uv.Y());
        if (at.SquareDistance(pole) < least) {
            least = at.SquareDistance(pole);
            nearest = u;
            foot = at;
        }
    }
    gp_Pnt edge;
    piece.rail.curve.D0(piece.rail.t(nearest), edge);
    return gp_Vec(foot, pole).Dot(gp_Vec(foot, edge)) > 0;
}

std::pair<Handle(Geom_Surface), TopLoc_Location> surfaceOf(const TopoDS_Face& face) {
    TopLoc_Location location;
    const auto surface = BRep_Tool::Surface(face, location);
    return {surface, location};
}

gp_Pnt2d nearPeriod(gp_Pnt2d p, const gp_Pnt2d& target, const Handle(Geom_Surface)& surface) {
    if (surface->IsUPeriodic()) p.SetX(loft::nearPeriod(p.X(), target.X(), surface->UPeriod()));
    if (surface->IsVPeriodic()) p.SetY(loft::nearPeriod(p.Y(), target.Y(), surface->VPeriod()));
    return p;
}

Handle(Geom2d_Curve) aligned(const Handle(Geom2d_Curve)& curve, double t, const gp_Pnt2d& expected,
                             const Handle(Geom_Surface)& surface) {
    const gp_Pnt2d at = curve->Value(t);
    const gp_Vec2d shift(at, nearPeriod(at, expected, surface));
    if (shift.SquareMagnitude() == 0) return curve;
    return Handle(Geom2d_Curve)::DownCast(curve->Translated(shift));
}

Handle(Geom2d_Curve) between(const Handle(Geom2d_Curve)& curve, double from, double to) { return new Geom2d_TrimmedCurve(curve, from, to); }

std::vector<Meeting> meetings(const Handle(Geom2d_Curve)& one, const Handle(Geom2d_Curve)& two, const Handle(Geom_Surface)& surface) {
    const double tu = surface->IsUPeriodic() ? surface->UPeriod() : 0, tv = surface->IsVPeriodic() ? surface->VPeriod() : 0;
    // The parameter tolerance relative to the parameters' size: where they run to 1e5 (a plane's
    // millimetres on a large body) a bare 1e-12 is below their rounding, and a meeting exactly at a
    // knot where the curve is only C0 falls between its pieces.
    double size = 1;
    for (const auto& curve : {one, two})
        for (const double t : {curve->FirstParameter(), curve->LastParameter()})
            if (!Precision::IsInfinite(t)) size = std::max({size, std::abs(curve->Value(t).X()), std::abs(curve->Value(t).Y())});
    const double tolerance = tol::parameter * size;
    std::vector<Meeting> out;
    for (int i = tu > 0 ? -1 : 0; i <= (tu > 0 ? 1 : 0); ++i)
        for (int j = tv > 0 ? -1 : 0; j <= (tv > 0 ? 1 : 0); ++j) {
            const gp_Vec2d shift(i * tu, j * tv);
            const auto laid = shift.SquareMagnitude() == 0 ? two : Handle(Geom2d_Curve)::DownCast(two->Translated(shift));
            Geom2dAPI_InterCurveCurve meet(one, laid, tolerance);
            for (int k = 1; k <= meet.NbPoints(); ++k) {
                const auto& point = meet.Intersector().Point(k);
                out.push_back({point.ParamOnFirst(), point.ParamOnSecond(), shift});
            }
        }
    return out;
}

std::optional<Meeting> nearest(const std::vector<Meeting>& found, double near) {
    std::optional<Meeting> best;
    for (const auto& m : found)
        if (!best || std::abs(m.first - near) < std::abs(best->first - near)) best = m;
    return best;
}

std::optional<OnCurve> onCurve(const gp_Pnt& p, const Handle(Geom_Curve)& curve) {
    std::optional<OnCurve> best;
    GeomAPI_ProjectPointOnCurve project(p, curve);
    if (project.NbPoints() > 0) best = OnCurve{project.LowerDistanceParameter(), project.LowerDistance()};
    for (const double t : {curve->FirstParameter(), curve->LastParameter()})
        if (!Precision::IsInfinite(t) && (!best || curve->Value(t).Distance(p) < best->distance)) best = OnCurve{t, curve->Value(t).Distance(p)};
    return best;
}

TopoDS_Edge edgeBetween(const Handle(Geom_Curve)& curve, const TopoDS_Vertex& first, double t0, const TopoDS_Vertex& second, double t1,
                        const std::vector<OnSurface>& on, double tolerance) {
    // A null vertex would be written through: in WASM that overwrites memory rather than traps.
    if (first.IsNull() || second.IsNull()) throw Failure("filletV2", no_solution, "a new edge has no vertex at an end", {});
    BRep_Builder builder;
    TopoDS_Edge edge;
    builder.MakeEdge(edge, curve, tolerance);
    for (const auto& c : on) builder.UpdateEdge(edge, c.curve, c.surface, c.location, tolerance);
    const bool rising = t0 < t1;
    const TopoDS_Vertex start = TopoDS::Vertex((rising ? first : second).Oriented(TopAbs_FORWARD));
    const TopoDS_Vertex end = TopoDS::Vertex((rising ? second : first).Oriented(TopAbs_REVERSED));
    builder.Add(edge, start);
    builder.Add(edge, end);
    builder.Range(edge, std::min(t0, t1), std::max(t0, t1));
    // Each end's parameter is the range's (BRep_Tool::Parameter takes it from the vertex's place in
    // the edge): a vertex is told only where it has no place on a curve yet, one just made, so that
    // the body's own vertices, shared with the input, are never changed.
    const auto fresh = [](const TopoDS_Vertex& v) {
        const auto tv = Handle(BRep_TVertex)::DownCast(v.TShape());
        return !tv.IsNull() && tv->Points().IsEmpty();
    };
    if (fresh(start)) builder.UpdateVertex(start, std::min(t0, t1), edge, tolerance);
    if (fresh(end)) builder.UpdateVertex(end, std::max(t0, t1), edge, tolerance);
    return edge;
}

Handle(Geom2d_Curve) laidOn(const Handle(Geom_Curve)& curve, double lo, double hi, const Handle(Geom_Surface)& surface) {
    const Handle(Geom2d_Curve) projected = GeomProjLib::Curve2d(curve, lo, hi, surface);
    if (projected.IsNull()) return projected;
    const auto off = [&](const Handle(Geom2d_Curve)& c) {
        double worst = 0;
        for (int i = 0; i <= 16; ++i) {
            const double t = lo + (hi - lo) * i / 16;
            const gp_Pnt2d p = c->Value(t);
            worst = std::max(worst, surface->Value(p.X(), p.Y()).Distance(curve->Value(t)));
        }
        return worst;
    };
    const double missed = off(projected);
    if (missed <= tol::exactLength) return projected;
    gp_Pnt p;
    gp_Vec du, dv;
    const gp_Pnt2d middle = projected->Value((lo + hi) / 2);
    surface->D1(middle.X(), middle.Y(), p, du, dv);
    const double speed = std::max(du.Magnitude(), dv.Magnitude());
    if (!(speed > limit::degenerate)) return projected;
    // Each point's own parameters, by Newton's method on the surface from the given curve's (within
    // about the length tolerance of them, quadratic convergence reaches rounding in three steps).
    bool lost = false;
    const auto own = [&](double t) {
        const gp_Pnt target = curve->Value(t);
        gp_Pnt2d uv = projected->Value(t);
        for (int step = 0; step < 3; ++step) {
            gp_Pnt at;
            gp_Vec su, sv;
            surface->D1(uv.X(), uv.Y(), at, su, sv);
            const gp_Vec off(at, target);
            const double a = su.Dot(su), b = su.Dot(sv), c = sv.Dot(sv), det = a * c - b * b;
            if (!(det > limit::degenerate * a * c)) {
                lost = true; // a pole or an apex of the surface
                break;
            }
            uv.SetCoord(uv.X() + (off.Dot(su) * c - off.Dot(sv) * b) / det, uv.Y() + (off.Dot(sv) * a - off.Dot(su) * b) / det);
        }
        return uv;
    };
    // A curve that itself lies off the surface by about as much (a fitted meeting's) gains nothing.
    double apart = 0;
    for (int i = 0; i <= 4; ++i) {
        const double t = lo + (hi - lo) * i / 4;
        const gp_Pnt2d uv = own(t);
        apart = std::max(apart, surface->Value(uv.X(), uv.Y()).Distance(curve->Value(t)));
    }
    if (lost || !(apart < missed / 2)) return projected;
    const auto laid = interpolated(own, lo, hi, tol::exactLength / speed);
    return lost || laid.IsNull() || !(off(laid) < missed) ? projected : Handle(Geom2d_Curve)(laid);
}

Handle(Geom2d_Curve) segment2d(const gp_Pnt2d& a, const gp_Pnt2d& b, double t0, double t1) {
    NCollection_Array1<gp_Pnt2d> poles(1, 2);
    poles(1) = a;
    poles(2) = b;
    NCollection_Array1<double> knots(1, 2);
    knots(1) = t0;
    knots(2) = t1;
    NCollection_Array1<int> multiplicities(1, 2);
    multiplicities.Init(2);
    return new Geom2d_BSplineCurve(poles, knots, multiplicities, 1);
}

std::optional<Meeting> crossing(const Handle(Geom2d_Curve)& contact, double from, double to, double near, const TopoDS_Edge& edge,
                                const TopoDS_Face& support) {
    double f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(edge, support, f, l);
    if (pcurve.IsNull()) return {};
    const auto surface = BRep_Tool::Surface(support);
    if (const auto found = nearest(meetings(between(contact, from, to), between(pcurve, f, l), surface), near)) return found;
    // The intersector can miss a crossing on a short trim of a contact that it finds on the whole.
    std::vector<Meeting> within;
    for (const auto& m : meetings(contact, between(pcurve, f, l), surface))
        if (m.first >= from && m.first <= to) within.push_back(m);
    return nearest(within, near);
}

/// How far from `vertex` the contact between `from` and `to` lies on its support, at most: a
/// place where it crosses an edge of the support through the vertex lies no further than that,
/// so no reach need be guessed (a wedge the stop face makes with the support acutely, where the
/// contact crosses far from the vertex, asks for more than a blend's size).
static double contactReach(const Handle(Geom2d_Curve)& contact, double from, double to, const TopoDS_Face& support, const TopoDS_Vertex& vertex) {
    const auto surface = BRep_Tool::Surface(support);
    const gp_Pnt at = BRep_Tool::Pnt(vertex);
    constexpr int samples = 16;
    double reach = 0;
    for (int i = 0; i <= samples; ++i) {
        const gp_Pnt2d uv = contact->Value(from + (to - from) * i / samples);
        reach = std::max(reach, surface->Value(uv.X(), uv.Y()).Distance(at));
    }
    return reach;
}

std::optional<Meeting> beyond(const Handle(Geom2d_Curve)& contact, double from, double to, double near, const TopoDS_Edge& edge,
                              const TopoDS_Vertex& vertex, const TopoDS_Face& support, const TopoDS_Face& stop, double tolerance) {
    const double reach = contactReach(contact, from, to, support, vertex);
    const TopoDS_Edge forward = TopoDS::Edge(edge.Oriented(TopAbs_FORWARD));
    double f, l, f2, l2;
    const auto stored = BRep_Tool::Curve(forward, f, l);
    const auto storedOnSupport = BRep_Tool::CurveOnSurface(forward, support, f2, l2);
    if (stored.IsNull() || storedOnSupport.IsNull()) return {};
    const auto curve = basisOf(stored);
    const auto pcurve = basisOf(storedOnSupport);
    const bool line = curve->IsKind(STANDARD_TYPE(Geom_Line)) && pcurve->IsKind(STANDARD_TYPE(Geom2d_Line));
    const bool circle = curve->IsKind(STANDARD_TYPE(Geom_Circle)) &&
                        (pcurve->IsKind(STANDARD_TYPE(Geom2d_Line)) || pcurve->IsKind(STANDARD_TYPE(Geom2d_Circle)));
    if (!line && !circle) return {};
    const bool atFirst = TopExp::FirstVertex(forward).IsSame(vertex);
    gp_Pnt p;
    gp_Vec d;
    curve->D1(atFirst ? f : l, p, d);
    if (!(d.Magnitude() > 0)) return {};
    double run = reach / d.Magnitude();
    if (circle) run = std::min(run, 2 * M_PI - (l - f));
    const double e0 = atFirst ? f - run : l, e1 = atFirst ? f : l + run;
    const auto best = nearest(meetings(between(contact, from, to), between(pcurve, e0, e1), BRep_Tool::Surface(support)), near);
    if (!best) return {};
    const gp_Pnt at = curve->Value(best->second);
    for (const auto& face : {support, stop}) {
        GeomAPI_ProjectPointOnSurf on(at, BRep_Tool::Surface(face));
        if (on.NbPoints() == 0 || on.LowerDistance() > tolerance) return {};
    }
    return best;
}

std::optional<Continuation> continuation(const Handle(Geom2d_Curve)& contact, double from, double to, double near,
                                         const TopoDS_Edge& edge, const TopoDS_Vertex& vertex, const TopoDS_Face& support,
                                         const TopoDS_Face& stop, double tolerance) {
    const double reach = contactReach(contact, from, to, support, vertex);
    const auto supportSurface = BRep_Tool::Surface(support), stopSurface = BRep_Tool::Surface(stop);
    const SurfaceMeeting meet(supportSurface, stopSurface, tolerance);
    const gp_Pnt at = BRep_Tool::Pnt(vertex);
    double f, l;
    const auto own = BRep_Tool::Curve(edge, f, l);
    if (own.IsNull()) return {};
    const gp_Pnt inside = own->Value((f + l) / 2);
    std::optional<Continuation> best;
    for (int i = 1; i <= meet.count(); ++i) {
        const auto line = meet.curve(i);
        const auto onLine = onCurve(at, line);
        if (!onLine || onLine->distance > tol::resolved(limit::together, reach, at.XYZ().Modulus())) continue;
        const double tVertex = onLine->parameter;
        const auto toEdge = onCurve(inside, line);
        if (!toEdge) continue;
        const double tEdge = toEdge->parameter;
        // Over the span that the contact can reach.
        gp_Pnt p;
        gp_Vec d;
        line->D1(tVertex, p, d);
        if (!(d.Magnitude() > 0)) continue;
        const double span = 2 * reach / d.Magnitude();
        const double lo = std::max(line->FirstParameter(), tVertex - span), hi = std::min(line->LastParameter(), tVertex + span);
        const auto onS1 = meet.on(i, 1, lo, hi), onS2 = meet.on(i, 2, lo, hi);
        if (onS1.IsNull() || onS2.IsNull()) continue;
        // Along the meeting from the vertex, the shorter way round where it closes (an ellipse
        // round a hole): a place a period off is the same place.
        const auto along = [&](double t) { return line->IsPeriodic() ? std::remainder(t - tVertex, line->Period()) : t - tVertex; };
        for (const auto& m : meetings(between(contact, from, to), onS1, supportSurface)) {
            // Past the vertex, away from the edge, and within reach of it.
            if (!(along(m.second) * along(tEdge) < 0) || line->Value(m.second).Distance(at) > reach) continue;
            if (best && std::abs(m.first - near) >= std::abs(best->u - near)) continue;
            const auto onSupport = m.shift.SquareMagnitude() == 0 ? onS1 : Handle(Geom2d_Curve)::DownCast(onS1->Translated(m.shift));
            best = Continuation{m.first, m.second, tVertex, line, onSupport, onS2};
        }
    }
    return best;
}

TopoDS_Shape withSolid(const TopoDS_Shape& body, const TopoDS_Shape& old, const TopoDS_Shape& solid) {
    if (body.IsSame(old)) return solid;
    if (body.ShapeType() != TopAbs_COMPOUND && body.ShapeType() != TopAbs_COMPSOLID) return body;
    BRep_Builder builder;
    TopoDS_Compound out;
    builder.MakeCompound(out);
    for (TopoDS_Iterator child(body); child.More(); child.Next()) builder.Add(out, withSolid(child.Value(), old, solid));
    return out;
}

Edit& editOf(Edits& edits, const TopoDS_Shape& face) {
    for (auto& edit : edits)
        if (edit.face.IsSame(face)) return edit;
    edits.push_back({TopoDS::Face(face), {}, {}});
    return edits.back();
}

Placed placedSheet(const Sheet& sheet) { return {sheet.surface, sheet.surface, {}}; }

Placed placedFace(const TopoDS_Face& face) {
    const auto [basis, location] = surfaceOf(face);
    return {BRep_Tool::Surface(face), basis, location};
}

/// Where the meeting comes in two branches, one through each corner, that end at one place (the
/// surfaces' meeting split where one of them is only G1, as a full round's two arcs are): one edge
/// along both, interpolated at places on them with the same parameter on its curve and on each
/// surface, and held to the branches between the places.
static std::optional<TopoDS_Edge> joinedMeeting(const SurfaceMeeting& meet, const Placed& one, const gp_Pnt2d& firstOnOne, const Placed& two,
                                                const gp_Pnt2d& firstOnTwo, const TopoDS_Vertex& first, const TopoDS_Vertex& second,
                                                double scale) {
    const gp_Pnt a = BRep_Tool::Pnt(first), b = BRep_Tool::Pnt(second);
    // OpenCascade's intersector is asked for the meeting to the length tolerance: a share of a
    // blend smaller than that cannot be told.
    const double together = std::max(tol::length, tol::resolved(limit::together, scale, a.XYZ().Modulus()));
    for (int i = 1; i <= meet.count(); ++i)
        for (int j = 1; j <= meet.count(); ++j) {
            if (i == j) continue;
            const auto ta = onCurve(a, meet.curve(i)), tb = onCurve(b, meet.curve(j));
            if (!ta || !tb || ta->distance > together || tb->distance > together) continue;
            // The ends of the two branches that meet.
            const auto& ci = meet.curve(i);
            const auto& cj = meet.curve(j);
            std::optional<std::pair<double, double>> joint;
            for (const double ei : {ci->FirstParameter(), ci->LastParameter()})
                for (const double ej : {cj->FirstParameter(), cj->LastParameter()})
                    if (ci->Value(ei).Distance(cj->Value(ej)) < together) joint = {ei, ej};
            if (!joint) continue;
            const auto oneI = meet.on(i, 1, std::min(ta->parameter, joint->first), std::max(ta->parameter, joint->first));
            const auto twoI = meet.on(i, 2, std::min(ta->parameter, joint->first), std::max(ta->parameter, joint->first));
            const auto oneJ = meet.on(j, 1, std::min(tb->parameter, joint->second), std::max(tb->parameter, joint->second));
            const auto twoJ = meet.on(j, 2, std::min(tb->parameter, joint->second), std::max(tb->parameter, joint->second));
            if (oneI.IsNull() || twoI.IsNull() || oneJ.IsNull() || twoJ.IsNull()) return {};
            // Each branch interpolated at places with the same parameter on its curve and on each
            // surface (each parameter curve taken by whole periods to pass where the first corner
            // is), over [0, 1/2] and [1/2, 1], and the two joined end to end: the parameter
            // curves' speeds may jump where the branches meet.
            const auto alignedOne = aligned(oneI, ta->parameter, firstOnOne, one.surface), alignedTwo = aligned(twoI, ta->parameter, firstOnTwo, two.surface);
            const Handle(Geom2d_Curve) branches[2][2] = {
                {alignedOne, alignedTwo},
                {aligned(oneJ, joint->second, alignedOne->Value(joint->first), one.surface), aligned(twoJ, joint->second, alignedTwo->Value(joint->first), two.surface)}};
            const Handle(Geom_Curve) lines[2] = {ci, cj};
            const double ends[2][2] = {{ta->parameter, joint->first}, {joint->second, tb->parameter}};
            // Places doubled until the curves between them lie on both surfaces.
            for (int count = 17; count <= 1025; count = 2 * count - 1) {
                std::vector<Handle(Geom_BSplineCurve)> curves;
                std::vector<Handle(Geom2d_BSplineCurve)> onSurface[2];
                double error = 0;
                for (int h = 0; h < 2; ++h) {
                    NCollection_Array1<gp_Pnt> points(1, count);
                    NCollection_Array1<gp_Pnt2d> uv[2] = {NCollection_Array1<gp_Pnt2d>(1, count), NCollection_Array1<gp_Pnt2d>(1, count)};
                    Handle(NCollection_HArray1<double>) parameters = new NCollection_HArray1<double>(1, count);
                    const auto at = [&](double share) { return ends[h][0] + (ends[h][1] - ends[h][0]) * share; };
                    for (int k = 0; k < count; ++k) {
                        const double t = at(double(k) / (count - 1));
                        gp_Pnt p;
                        gp_Vec d;
                        lines[h]->D1(t, p, d);
                        gp_Pnt2d q[2] = {branches[h][0]->Value(t), branches[h][1]->Value(t)};
                        // The branch's place taken onto both surfaces exactly (OpenCascade's meeting is an
                        // approximation): where they meet in the plane across the branch there.
                        std::array<double, 4> x = {q[0].X(), q[0].Y(), q[1].X(), q[1].Y()};
                        const auto residual = [&](const std::array<double, 4>& y, std::array<std::array<double, 5>, 4>& m, double& norm) {
                            gp_Pnt p0, p1;
                            gp_Vec du0, dv0, du1, dv1;
                            one.surface->D1(y[0], y[1], p0, du0, dv0);
                            two.surface->D1(y[2], y[3], p1, du1, dv1);
                            const gp_Vec f(p1, p0);
                            for (int i = 0; i < 3; ++i) m[i] = {du0.Coord(i + 1), dv0.Coord(i + 1), -du1.Coord(i + 1), -dv1.Coord(i + 1), -f.Coord(i + 1)};
                            const double g = gp_Vec(p, p0).Dot(d);
                            m[3] = {du0.Dot(d), dv0.Dot(d), 0, 0, -g};
                            norm = f.Magnitude() + std::abs(g) / std::max(d.Magnitude(), limit::degenerate);
                            return std::isfinite(norm);
                        };
                        if (newton(x, tol::resolved(limit::solved, scale, p.XYZ().Modulus()), residual)) {
                            q[0].SetCoord(x[0], x[1]);
                            q[1].SetCoord(x[2], x[3]);
                            p = one.surface->Value(x[0], x[1]);
                        }
                        points(k + 1) = p;
                        for (int surface = 0; surface < 2; ++surface) uv[surface](k + 1) = q[surface];
                        parameters->SetValue(k + 1, (h + double(k) / (count - 1)) / 2);
                    }
                    if (h == 0) points(1) = a;
                    else points(count) = b;
                    GeomAPI_Interpolate curve(new NCollection_HArray1<gp_Pnt>(points), parameters, false, tol::parameter);
                    curve.Perform();
                    if (!curve.IsDone()) return {};
                    curves.push_back(curve.Curve());
                    for (int surface = 0; surface < 2; ++surface) {
                        Geom2dAPI_Interpolate fit(new NCollection_HArray1<gp_Pnt2d>(uv[surface]), parameters, false, tol::parameter);
                        fit.Perform();
                        if (!fit.IsDone()) return {};
                        onSurface[surface].push_back(fit.Curve());
                    }
                    // Between the places, against both surfaces.
                    const Placed* placed[2] = {&one, &two};
                    for (int k = 0; k + 1 < count; ++k) {
                        const double share = (k + 0.5) / (count - 1), u = (h + share) / 2;
                        const gp_Pnt p = curves.back()->Value(u);
                        for (int surface = 0; surface < 2; ++surface) {
                            const gp_Pnt2d q = onSurface[surface].back()->Value(u);
                            error = std::max(error, placed[surface]->surface->Value(q.X(), q.Y()).Distance(p));
                        }
                    }
                }
                if (error > together) continue;
                const auto curve = endToEnd(curves, together);
                const auto curveOne = endToEnd(onSurface[0], tol::parameterMerge), curveTwo = endToEnd(onSurface[1], tol::parameterMerge);
                if (curve.IsNull() || curveOne.IsNull() || curveTwo.IsNull()) return {};
                return edgeBetween(curve, first, 0, second, 1, {{curveOne, one.basis, one.location}, {curveTwo, two.basis, two.location}},
                                   std::max(tol::length, error));
            }
            return {};
        }
    return {};
}

/// Where an extruded sheet whose section is polynomial meets a plane across the extrusion: the
/// section carried along the extrusion onto the plane, which is an affine image of it (its poles
/// moved, the same degree and knots), with its curve on the sheet (the section's own parameter,
/// and the distance along the extrusion as a spline of the same knots) and on the plane, all exact.
/// GeomInt approximates this meeting by a C1 curve of degree 8, whose knots a later blend along
/// the edge takes for breaks of its rail (and is only C0 there). Empty where the surfaces are not
/// these.
static std::optional<TopoDS_Edge> extrusionMeetsPlane(const Placed& one, const Placed& two, const TopoDS_Vertex& first,
                                                      const TopoDS_Vertex& second, double scale) {
    for (const auto& [sheet, face] : {std::pair{one, two}, std::pair{two, one}}) {
        const auto extrusion = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(sheet.surface);
        const auto plane = Handle(Geom_Plane)::DownCast(face.surface);
        if (extrusion.IsNull() || plane.IsNull()) continue;
        const auto section = Handle(Geom_BSplineCurve)::DownCast(extrusion->BasisCurve());
        if (section.IsNull() || section->IsRational()) return {};
        const gp_Vec along(extrusion->Direction());
        const gp_Pln pln = plane->Pln();
        const gp_Vec normal(pln.Axis().Direction());
        // Across the plane, not along it.
        if (std::abs(along.Dot(normal)) <= tol::coarseAngle) return {};
        const int n = section->NbPoles();
        NCollection_Array1<gp_Pnt> poles(1, n);
        NCollection_Array1<gp_Pnt2d> onSheet(1, n);
        NCollection_Array1<double> flat(1, section->NbPoles() + section->Degree() + 1);
        section->KnotSequence(flat);
        for (int i = 1; i <= n; ++i) {
            const gp_Pnt p = section->Pole(i);
            const double v = gp_Vec(p, pln.Location()).Dot(normal) / along.Dot(normal);
            poles(i) = p.Translated(along * v);
            // The section's parameter is the sum of its basis functions times their Greville
            // abscissae.
            double greville = 0;
            for (int j = 1; j <= section->Degree(); ++j) greville += flat(i + j);
            onSheet(i) = gp_Pnt2d(greville / section->Degree(), v);
        }
        NCollection_Array1<double> knots(1, section->NbKnots());
        NCollection_Array1<int> multiplicities(1, section->NbKnots());
        section->Knots(knots);
        section->Multiplicities(multiplicities);
        const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, section->Degree());
        const Handle(Geom2d_Curve) onSheetCurve = new Geom2d_BSplineCurve(onSheet, knots, multiplicities, section->Degree());
        const Handle(Geom2d_Curve) onPlane = GeomAPI::To2d(curve, pln);
        if (onPlane.IsNull()) return {};
        GeomAPI_ProjectPointOnCurve a(BRep_Tool::Pnt(first), curve), b(BRep_Tool::Pnt(second), curve);
        if (a.NbPoints() == 0 || b.NbPoints() == 0) return {};
        const double together = limit::together * scale;
        if (a.LowerDistance() > together || b.LowerDistance() > together) return {};
        return edgeBetween(curve, first, a.LowerDistanceParameter(), second, b.LowerDistanceParameter(),
                           {{onSheetCurve, sheet.basis, sheet.location}, {onPlane, face.basis, face.location}},
                           std::max({tol::length, a.LowerDistance(), b.LowerDistance()}));
    }
    return {};
}

static std::optional<TopoDS_Edge> crossedSurfaces(const Placed& one, const gp_Pnt2d& oneFirst, const gp_Pnt2d& oneSecond,
                                                  const Placed& two, const gp_Pnt2d& twoFirst, const gp_Pnt2d& twoSecond,
                                                  const TopoDS_Vertex& first, const TopoDS_Vertex& second);

std::optional<TopoDS_Edge> meetingEdge(const Placed& one, const gp_Pnt2d& firstOnOne, const Placed& two, const gp_Pnt2d& firstOnTwo,
                                       const TopoDS_Vertex& first, const TopoDS_Vertex& second, double scale) {
    if (auto exact = extrusionMeetsPlane(one, two, first, second, scale)) return exact;
    const double tolerance = tol::length;
    // A small blend's meeting laid on its surfaces as true to its size as a large one's (as the fit).
    const SurfaceMeeting meet(one.surface, two.surface, tolerance, tolerance * std::min(1.0, scale));
    const gp_Pnt a = BRep_Tool::Pnt(first), b = BRep_Tool::Pnt(second);
    int best = 0;
    double gap = HUGE_VAL, pa = 0, pb = 0;
    for (int i = 1; i <= meet.count(); ++i) {
        const auto ta = onCurve(a, meet.curve(i)), tb = onCurve(b, meet.curve(i));
        if (!ta || !tb) continue;
        const double g = std::max(ta->distance, tb->distance);
        if (g < gap) {
            gap = g;
            best = i;
            pa = ta->parameter;
            pb = tb->parameter;
        }
    }
    // The two known corners can seed the meeting when the general intersector does not return a
    // complete branch. The chord sections below certify the whole connecting edge.
    const auto seeded = [&]() -> std::optional<TopoDS_Edge> {
        gp_Pnt2d last[2];
        const Placed* surfaces[2] = {&one, &two};
        const gp_Pnt2d initial[2] = {firstOnOne, firstOnTwo};
        for (int i = 0; i < 2; ++i) {
            GeomAPI_ProjectPointOnSurf project(b, surfaces[i]->surface);
            if (project.NbPoints() == 0 || project.LowerDistance() > tolerance) return {};
            double u, v;
            project.LowerDistanceParameters(u, v);
            last[i] = nearPeriod(gp_Pnt2d(u, v), initial[i], surfaces[i]->surface);
            if (surfaces[i]->surface->Value(initial[i].X(), initial[i].Y()).Distance(a) > tolerance) return {};
        }
        return crossedSurfaces(one, firstOnOne, last[0], two, firstOnTwo, last[1], first, second);
    };
    if (best == 0 || gap > std::max(tol::length, tol::resolved(limit::together, scale, a.XYZ().Modulus()))) {
        if (auto joined = joinedMeeting(meet, one, firstOnOne, two, firstOnTwo, first, second, scale)) return joined;
        return seeded();
    }
    // A blend's end turns less than half way round. On a closed curve that is not periodic (two
    // cylinders' meeting, a B-spline loop) the shorter way may pass its ends: where one place is
    // at an end, it is taken at the other; through the loop's ends elsewhere it is found between its
    // corners instead.
    const auto& line = meet.curve(best);
    if (line->IsPeriodic() && std::abs(pb - pa) > line->Period() / 2) (pa < pb ? pa : pb) += line->Period();
    else if (!line->IsPeriodic() && line->IsClosed()) {
        const double f = line->FirstParameter(), l = line->LastParameter(), near = tol::parameterMerge * (l - f);
        if (std::abs(pb - pa) > (l - f) / 2) {
            double& atEnd = std::min(std::abs(pa - f), std::abs(pa - l)) <= near ? pa : pb;
            if (std::min(std::abs(atEnd - f), std::abs(atEnd - l)) > near) return seeded();
            atEnd = std::abs(atEnd - f) <= near ? l : f;
        }
    }
    const double lo = std::min(pa, pb), hi = std::max(pa, pb);
    const auto onOne = meet.on(best, 1, lo, hi), onTwo = meet.on(best, 2, lo, hi);
    if (onOne.IsNull() || onTwo.IsNull()) return {};
    return edgeBetween(line, first, pa, second, pb,
                       {{aligned(onOne, pa, firstOnOne, one.surface), one.basis, one.location},
                        {aligned(onTwo, pa, firstOnTwo, two.surface), two.basis, two.location}},
                       std::max(tolerance, gap));
}

namespace {
/// A place where two sheets meet: on the first sheet's surface, on the second's, and in space.
struct Place {
    gp_Pnt2d a, b;
    gp_Pnt p;
};

/// The edge through the places where two sheets meet from `first` to `last`, `solve(t, place)`
/// finding the place at t in (0, 1) from the place before it, moved on as the straight line between
/// the ends would: the curve and both parameter curves interpolate the same places at the same
/// parameter, doubled until the places between them are met within the length tolerance.
std::optional<TopoDS_Edge> throughPlaces(const Placed& one, const Placed& two, const Place& first,
                                         const Place& last, const std::function<bool(double, Place&)>& solve, const TopoDS_Vertex& from,
                                         const TopoDS_Vertex& to) {
    const double tolerance = tol::length;
    for (int count = 16; count <= 1024; count *= 2) {
        std::vector<Place> places(count + 1);
        places.front() = first;
        places.back() = last;
        bool solved = true;
        for (int i = 1; i < count && solved; ++i) {
            Place& place = places[i];
            place = places[i - 1];
            const double f = 1.0 / (count - i + 1);
            place.a = gp_Pnt2d(place.a.XY() + (last.a.XY() - place.a.XY()) * f);
            place.b = gp_Pnt2d(place.b.XY() + (last.b.XY() - place.b.XY()) * f);
            solved = solve(double(i) / count, place);
        }
        if (!solved) continue;
        // The last solved place must run on to the second corner: they may touch there.
        NCollection_Array1<gp_Pnt> points(1, count + 1);
        NCollection_Array1<gp_Pnt2d> onOne(1, count + 1), onTwo(1, count + 1);
        Handle(NCollection_HArray1<double>) parameters = new NCollection_HArray1<double>(1, count + 1);
        for (int i = 0; i <= count; ++i) {
            points(i + 1) = places[i].p;
            onOne(i + 1) = places[i].a;
            onTwo(i + 1) = places[i].b;
            parameters->SetValue(i + 1, double(i) / count);
        }
        const auto interpolate3 = [&](const NCollection_Array1<gp_Pnt>& at) {
            Handle(NCollection_HArray1<gp_Pnt>) values = new NCollection_HArray1<gp_Pnt>(at);
            GeomAPI_Interpolate fit(values, parameters, false, tol::parameter);
            fit.Perform();
            return fit.IsDone() ? Handle(Geom_Curve)(fit.Curve()) : Handle(Geom_Curve)();
        };
        const auto interpolate2 = [&](const NCollection_Array1<gp_Pnt2d>& at) {
            Handle(NCollection_HArray1<gp_Pnt2d>) values = new NCollection_HArray1<gp_Pnt2d>(at);
            Geom2dAPI_Interpolate fit(values, parameters, false, tol::parameter);
            fit.Perform();
            return fit.IsDone() ? Handle(Geom2d_Curve)(fit.Curve()) : Handle(Geom2d_Curve)();
        };
        const auto curve = interpolate3(points);
        const auto curveOne = interpolate2(onOne), curveTwo = interpolate2(onTwo);
        if (curve.IsNull() || curveOne.IsNull() || curveTwo.IsNull()) return {};
        // Between the places, against the exact meeting there.
        double error = 0;
        for (int i = 0; i < count && error <= tolerance; ++i) {
            Place mid = places[i];
            const double t = (i + 0.5) / count;
            if (!solve(t, mid)) {
                error = HUGE_VAL;
                break;
            }
            const gp_Pnt2d a = curveOne->Value(t), b = curveTwo->Value(t);
            error = std::max({error, curve->Value(t).Distance(mid.p), one.surface->Value(a.X(), a.Y()).Distance(mid.p), two.surface->Value(b.X(), b.Y()).Distance(mid.p)});
        }
        if (error > tolerance) continue;
        return edgeBetween(curve, from, 0, to, 1, {{curveOne, one.basis, one.location}, {curveTwo, two.basis, two.location}}, tolerance);
    }
    return {};
}
} // namespace

std::optional<TopoDS_Edge> marchedMeeting(const Sheet& one, const gp_Pnt2d& oneFirst, const gp_Pnt2d& oneSecond, const Sheet& two,
                                          const gp_Pnt2d& twoFirst, const gp_Pnt2d& twoSecond, const TopoDS_Vertex& first,
                                          const TopoDS_Vertex& second) {
    const double tolerance = tol::length;
    // Marched along its length (Patrikalakis & Maekawa): each place is where the sheets meet on
    // the plane across the meeting a step on from the last, along the line where their tangent
    // planes meet, all four parameters solved by Newton's method. A step of length stays regular
    // where a sheet's own parameters do not: a mitre's meeting turns tightly by a narrow notch
    // between its walls, and it runs nearly along a contact where the contacts cross at a shallow
    // angle (a sheet's across parameter turns back there).
    // On the plane through `origin` across `normal`, from `place`; within the length it goes on
    // while it still gains, at most two steps more (to rounding): where the sheets nearly touch,
    // a residual only within the length leaves a place that over the angle between them off the
    // meeting.
    const auto solve = [&](const gp_Pnt& origin, const gp_Vec& normal, Place& place) {
        double x[4] = {place.a.X(), place.a.Y(), place.b.X(), place.b.Y()};
        std::optional<Place> found;
        double best = std::max(limit::marched, tol::rounding(origin.XYZ().Modulus()));
        int polished = 0;
        for (int iteration = 0; iteration < 50; ++iteration) {
            gp_Pnt pa, pb;
            gp_Vec dau, dav, dbu, dbv;
            one.surface->D1(x[0], x[1], pa, dau, dav);
            two.surface->D1(x[2], x[3], pb, dbu, dbv);
            const gp_Vec f(pb, pa);
            const double off = gp_Vec(origin, pa).Dot(normal);
            const double residual = std::max(f.Magnitude(), std::abs(off));
            if (residual < best) {
                found = Place{gp_Pnt2d(x[0], x[1]), gp_Pnt2d(x[2], x[3]), pa};
                best = residual;
                if (polished++ == 2) break;
            } else if (found)
                break;
            math_Matrix m(1, 4, 1, 4);
            math_Vector rhs(1, 4), step(1, 4);
            for (int i = 1; i <= 3; ++i) {
                m(i, 1) = dau.Coord(i);
                m(i, 2) = dav.Coord(i);
                m(i, 3) = -dbu.Coord(i);
                m(i, 4) = -dbv.Coord(i);
                rhs(i) = -f.Coord(i);
            }
            m(4, 1) = dau.Dot(normal);
            m(4, 2) = dav.Dot(normal);
            m(4, 3) = m(4, 4) = 0;
            rhs(4) = -off;
            math_Gauss gauss(m);
            if (!gauss.IsDone()) break;
            gauss.Solve(rhs, step);
            for (int i = 0; i < 4; ++i) x[i] += step(i + 1);
        }
        if (!found) return false;
        place = *found;
        return true;
    };
    // The meeting's direction at a place, the way `towards` goes.
    const auto along = [&](const Place& place, const gp_Vec& towards, gp_Vec& out) {
        gp_Pnt p;
        gp_Vec dau, dav, dbu, dbv;
        one.surface->D1(place.a.X(), place.a.Y(), p, dau, dav);
        two.surface->D1(place.b.X(), place.b.Y(), p, dbu, dbv);
        const gp_Vec d = dau.Crossed(dav).Crossed(dbu.Crossed(dbv));
        if (!(d.Magnitude() > limit::degenerate * dau.Magnitude() * dav.Magnitude() * dbu.Magnitude() * dbv.Magnitude())) return false;
        out = d.Normalized() * (d.Dot(towards) < 0 ? -1.0 : 1.0);
        return true;
    };
    const gp_Pnt A = BRep_Tool::Pnt(first), B = BRep_Tool::Pnt(second);
    double length = A.Distance(B);
    if (!(length > limit::degenerate)) return {};
    // Halving the steps takes a smooth meeting's error down some 16 times (cubic interpolation);
    // as the fit's refinement (blend2_fit.cpp), an error still a thousand tolerances that does
    // not halve in three is a break in it, and the meeting is declined; near the tolerance a slow
    // tail goes on. The places are bounded by memory (a long meeting at a small tolerance, as a
    // large body has, takes many).
    std::vector<double> errors;
    const auto stalled = [&](double error) {
        errors.push_back(error);
        return errors.size() > 3 && error > 1e3 * tolerance && !(error < errors[errors.size() - 4] / 2);
    };
    for (int count = 16; count <= (1 << 16); count *= 2) {
        const double h = length / count;
        std::vector<Place> places = {{oneFirst, twoFirst, A}};
        gp_Vec direction(A, B);
        bool solved = true;
        while (solved && places.back().p.Distance(B) > 1.5 * h) {
            // A step on along the meeting, from the last place moved on as the last step went (or
            // as the straight line between the corners would, at the first).
            if (places.size() > 4 * size_t(count) || !along(places.back(), places.size() == 1 ? gp_Vec(A, B) : direction, direction)) {
                solved = false;
                break;
            }
            Place place = places.back();
            if (places.size() > 1) {
                const Place& before = places[places.size() - 2];
                place.a = gp_Pnt2d(2 * place.a.XY() - before.a.XY());
                place.b = gp_Pnt2d(2 * place.b.XY() - before.b.XY());
            } else {
                place.a = gp_Pnt2d(place.a.XY() + (oneSecond.XY() - place.a.XY()) / count);
                place.b = gp_Pnt2d(place.b.XY() + (twoSecond.XY() - place.b.XY()) / count);
            }
            solved = solve(places.back().p.Translated(direction * h), direction, place);
            if (solved) places.push_back(place);
        }
        if (!solved) {
            if (stalled(HUGE_VAL)) break;
            continue;
        }
        // The last solved place runs on to the second corner: they may touch there. The places
        // are interpolated at their chord lengths.
        places.push_back({oneSecond, twoSecond, B});
        const int n = int(places.size());
        NCollection_Array1<gp_Pnt> points(1, n);
        NCollection_Array1<gp_Pnt2d> onOne(1, n), onTwo(1, n);
        Handle(NCollection_HArray1<double>) parameters = new NCollection_HArray1<double>(1, n);
        std::vector<double> chord(n, 0.0);
        for (int i = 1; i < n; ++i) chord[i] = chord[i - 1] + places[i].p.Distance(places[i - 1].p);
        for (int i = 0; i < n; ++i) {
            points(i + 1) = places[i].p;
            onOne(i + 1) = places[i].a;
            onTwo(i + 1) = places[i].b;
            parameters->SetValue(i + 1, i == n - 1 ? 1.0 : chord[i] / chord[n - 1]);
        }
        length = chord[n - 1];
        const auto interpolate3 = [&](const NCollection_Array1<gp_Pnt>& at) {
            Handle(NCollection_HArray1<gp_Pnt>) values = new NCollection_HArray1<gp_Pnt>(at);
            GeomAPI_Interpolate fit(values, parameters, false, tol::parameter);
            fit.Perform();
            return fit.IsDone() ? Handle(Geom_Curve)(fit.Curve()) : Handle(Geom_Curve)();
        };
        const auto interpolate2 = [&](const NCollection_Array1<gp_Pnt2d>& at) {
            Handle(NCollection_HArray1<gp_Pnt2d>) values = new NCollection_HArray1<gp_Pnt2d>(at);
            Geom2dAPI_Interpolate fit(values, parameters, false, tol::parameter);
            fit.Perform();
            return fit.IsDone() ? Handle(Geom2d_Curve)(fit.Curve()) : Handle(Geom2d_Curve)();
        };
        const auto curve = interpolate3(points);
        const auto curveOne = interpolate2(onOne), curveTwo = interpolate2(onTwo);
        if (curve.IsNull() || curveOne.IsNull() || curveTwo.IsNull()) return {};
        // Between the places, against the exact meeting there (across the chord's middle): the
        // curve's nearest place to it, and both parameter curves there on the curve.
        double error = 0;
        for (int i = 1; i < n - 2 && error <= tolerance; ++i) {
            Place mid = places[i];
            mid.a = gp_Pnt2d((places[i].a.XY() + places[i + 1].a.XY()) / 2);
            mid.b = gp_Pnt2d((places[i].b.XY() + places[i + 1].b.XY()) / 2);
            const gp_Vec span(places[i].p, places[i + 1].p);
            if (!(span.Magnitude() > limit::degenerate) ||
                !solve(gp_Pnt((places[i].p.XYZ() + places[i + 1].p.XYZ()) / 2), span.Normalized(), mid)) {
                error = HUGE_VAL;
                break;
            }
            const double lo = parameters->Value(i + 1), hi = parameters->Value(i + 2);
            double w = (lo + hi) / 2;
            for (int k = 0; k < 4; ++k) {
                gp_Pnt c;
                gp_Vec d1, d2;
                curve->D2(w, c, d1, d2);
                const gp_Vec off(mid.p, c);
                const double slope = d1.SquareMagnitude() + off.Dot(d2);
                if (!(slope > 0)) break;
                w = std::clamp(w - off.Dot(d1) / slope, lo, hi);
            }
            // Off the meeting the curve opens a gap between the sheets of that distance times the
            // sine of their angle: where they touch (a mitre's inner end) the meeting itself is
            // only as certain as the sheets, and any place between them is on both.
            gp_Pnt p;
            gp_Vec dau, dav, dbu, dbv;
            one.surface->D1(mid.a.X(), mid.a.Y(), p, dau, dav);
            two.surface->D1(mid.b.X(), mid.b.Y(), p, dbu, dbv);
            const gp_Vec na = dau.Crossed(dav), nb = dbu.Crossed(dbv);
            const double sine = na.Crossed(nb).Magnitude() / std::max(na.Magnitude() * nb.Magnitude(), std::numeric_limits<double>::min());
            const gp_Pnt c = curve->Value(w);
            const gp_Pnt2d a = curveOne->Value(w), b = curveTwo->Value(w);
            error = std::max({error, c.Distance(mid.p) * sine, one.surface->Value(a.X(), a.Y()).Distance(c), two.surface->Value(b.X(), b.Y()).Distance(c)});
        }
        if (error <= tolerance)
            return edgeBetween(curve, first, 0, second, 1, {{curveOne, one.surface, {}}, {curveTwo, two.surface, {}}}, tolerance);
        if (stalled(error)) break;
    }
    return {};
}

static std::optional<TopoDS_Edge> crossedSurfaces(const Placed& one, const gp_Pnt2d& oneFirst, const gp_Pnt2d& oneSecond,
                                                  const Placed& two, const gp_Pnt2d& twoFirst, const gp_Pnt2d& twoSecond,
                                                  const TopoDS_Vertex& first, const TopoDS_Vertex& second) {
    // A place at `t`: where the plane across the chord at that share of it meets both sheets, all
    // four parameters solved together.
    const gp_Pnt P = BRep_Tool::Pnt(first), Q = BRep_Tool::Pnt(second);
    const double length = P.Distance(Q);
    if (!(length > tol::length)) return {};
    const gp_Vec chord = gp_Vec(P, Q) / length;
    const auto solve = [&](double t, Place& place) {
        std::array<double, 4> x = {place.a.X(), place.a.Y(), place.b.X(), place.b.Y()};
        const auto residual = [&](const std::array<double, 4>& y, std::array<std::array<double, 5>, 4>& m, double& norm) {
            gp_Pnt pa, pb;
            gp_Vec dau, dav, dbu, dbv;
            one.surface->D1(y[0], y[1], pa, dau, dav);
            two.surface->D1(y[2], y[3], pb, dbu, dbv);
            const gp_Vec f(pb, pa);
            const double along = gp_Vec(P, pa).Dot(chord) - t * length;
            for (int i = 0; i < 3; ++i) m[i] = {dau.Coord(i + 1), dav.Coord(i + 1), -dbu.Coord(i + 1), -dbv.Coord(i + 1), -f.Coord(i + 1)};
            m[3] = {dau.Dot(chord), dav.Dot(chord), 0, 0, -along};
            norm = std::sqrt(f.SquareMagnitude() + along * along);
            return true;
        };
        if (!newton<4>(x, std::max(limit::marched, tol::rounding(P.XYZ().Modulus())), residual)) return false;
        place.a = gp_Pnt2d(x[0], x[1]);
        place.b = gp_Pnt2d(x[2], x[3]);
        place.p = one.surface->Value(x[0], x[1]);
        return true;
    };
    return throughPlaces(one, two, {oneFirst, twoFirst, P}, {oneSecond, twoSecond, Q}, solve, first, second);
}

std::optional<TopoDS_Edge> crossedMeeting(const Sheet& one, const gp_Pnt2d& oneFirst, const gp_Pnt2d& oneSecond, const Sheet& two,
                                          const gp_Pnt2d& twoFirst, const gp_Pnt2d& twoSecond, const TopoDS_Vertex& first,
                                          const TopoDS_Vertex& second) {
    return crossedSurfaces(placedSheet(one), oneFirst, oneSecond, placedSheet(two), twoFirst, twoSecond, first, second);
}

TopoDS_Edge isoSection(const Sheet& sheet, double rail, const TopoDS_Vertex& first, const TopoDS_Vertex& second) {
    return edgeBetween(sheet.section(rail), first, sheet.across[0], second, sheet.across[1], {}, tol::length);
}

bool leaving(const TopoDS_Edge& edge, const TopoDS_Vertex& vertex, gp_Vec& out, int end) {
    if (BRep_Tool::Degenerated(edge)) return false;
    TopoDS_Vertex first, last;
    TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), first, last);
    const BRepAdaptor_Curve curve(edge);
    gp_Pnt p;
    if (end != 1 && first.IsSame(vertex)) curve.D1(curve.FirstParameter(), p, out);
    else if (end != 0 && last.IsSame(vertex)) {
        curve.D1(curve.LastParameter(), p, out);
        out.Reverse();
    } else
        return false;
    if (!(out.Magnitude() > limit::degenerate)) return false;
    out.Normalize();
    return true;
}

bool tangent(const gp_Vec& a, const gp_Vec& b) { return a.Dot(b) < 0 && a.Crossed(b).Magnitude() <= tol::coarseAngle; }

bool normalAt(const TopoDS_Face& face, const TopoDS_Edge& edge, const TopoDS_Vertex& vertex, gp_Vec& out) {
    double first, last;
    const auto pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
    if (pcurve.IsNull()) return false;
    const gp_Pnt2d uv = pcurve->Value(BRep_Tool::Parameter(vertex, edge, face));
    const Side side{face, BRep_Tool::Surface(face), face.Orientation() == TopAbs_REVERSED ? -1.0 : 1.0, 0};
    Jet j;
    if (!jetAt(side, uv.X(), uv.Y(), j)) return false;
    out = j.n;
    return true;
}

bool smoothAt(const TopoDS_Face& a, const TopoDS_Face& b, const TopoDS_Edge& edge, const TopoDS_Vertex& vertex) {
    gp_Vec na, nb;
    return !a.IsNull() && !b.IsNull() && normalAt(a, edge, vertex, na) && normalAt(b, edge, vertex, nb) && na.Dot(nb) > 0 &&
           na.Crossed(nb).Magnitude() <= tol::coarseAngle;
}

} // namespace loft::blend2::detail

