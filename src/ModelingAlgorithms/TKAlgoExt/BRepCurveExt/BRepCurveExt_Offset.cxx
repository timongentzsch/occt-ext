// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <GProp_GProps.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Map.hxx>
#include <Precision.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace BRepExt {
namespace {
// How the gap at a corner of an offset is closed: 0 by an arc about the corner (Round), 1 by
// the chord of that arc (Linear), 2 by carrying the two neighbours on to where they cross
// (Natural). OpenCascade has arcs for every curve and crossings for lines and circles only, so
// Linear, and Natural beside any other curve, are made from the arcs.
//
// Everything here works on a copy of the curve moved to the origin and scaled to this size, and
// the result is taken back. OpenCascade's offset works to lengths of its own (1e-7 for "the
// same point"): what it returns for one curve differs with the curve's scale (a spline blob
// trimmed rightly at sizes 8 and 80 crosses itself at 800 and is refused at 80000) and, further
// from the origin, with its place. At one size every length below is a share of the curve.
constexpr double workSize = 10;
// Two points this near are one (OpenCascade's own measure, at the working size).
constexpr double same = Precision::Confusion();
// Ends of consecutive pieces this near are joined by a vertex that covers the gap. OpenCascade
// leaves gaps of 1.2e-7 between the pieces of a spline's offset and holds them by the vertex.
constexpr double joinGap = 100 * same;
// OpenCascade draws the offset of a spline to about 1e-4 of the distance (0.29997 for 0.3 was
// measured); a piece ten times that nearer to the source than the distance is not at the
// distance.
constexpr double nearShare = 1e-3;
// Two directions this near (radians) are one: the arcs OpenCascade joins pieces by are tangent
// to them to about 1e-9.
constexpr double tangent = BRepExt::tolerance::coarseAngle;

// A stretch of a curve as the offset runs along it, from `from` to `to` (from > to: backwards).
struct Run {
    Handle(Geom_Curve) curve;
    double from, to;
    bool join;
};

gp_Pnt at(const Run &run, bool end) { return run.curve->Value(end ? run.to : run.from); }

// The direction of travel at an end of a run.
gp_Vec heading(const Run &run, bool end) {
    gp_Pnt p;
    gp_Vec d;
    run.curve->D1(end ? run.to : run.from, p, d);
    return run.from > run.to ? -d : d;
}

// A run as an edge with ends of its own.
TopoDS_Edge edgeOf(const char *op, const Run &run) {
    BRepBuilderAPI_MakeEdge edge(run.curve, std::min(run.from, run.to), std::max(run.from, run.to));
    if (!edge.IsDone())
        fail(op, no_solution, "a piece of the offset could not be built");
    return run.from > run.to ? TopoDS::Edge(edge.Edge().Reversed()) : edge.Edge();
}

std::vector<Run> runsOf(const char *op, const TopoDS_Wire &wire,
                        const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> &joins) {
    std::vector<Run> runs;
    for (BRepTools_WireExplorer it(wire); it.More(); it.Next()) {
        double first, last;
        const auto curve = BRep_Tool::Curve(it.Current(), first, last);
        if (curve.IsNull())
            fail(op, no_solution, "the offset has an edge without a curve");
        const bool back = it.Current().Orientation() == TopAbs_REVERSED;
        runs.push_back({curve, back ? last : first, back ? first : last, joins.Contains(it.Current())});
    }
    return runs;
}

// The runs as one wire. Consecutive runs share a vertex that covers what lies between their
// ends: built with ends of their own, pieces 1.2e-7 apart do not make a wire. `apart` is what
// to say where two are further apart than a vertex may cover.
TopoDS_Wire wireOf(const char *op, const std::vector<Run> &runs, bool closed, const char *apart) {
    const size_t n = runs.size();
    BRep_Builder builder;
    std::vector<TopoDS_Vertex> starts(n), ends(n);
    for (size_t j = 0; j < n; ++j) {
        if (j + 1 == n && !closed)
            break;
        const size_t k = (j + 1) % n;
        const gp_Pnt a = at(runs[j], true), b = at(runs[k], false);
        if (a.Distance(b) > joinGap)
            fail(op, no_solution, apart);
        TopoDS_Vertex vertex;
        builder.MakeVertex(vertex, gp_Pnt((a.XYZ() + b.XYZ()) / 2), a.Distance(b) / 2 + same);
        ends[j] = starts[k] = vertex;
    }
    if (!closed) {
        starts[0] = BRepBuilderAPI_MakeVertex(at(runs[0], false)).Vertex();
        ends[n - 1] = BRepBuilderAPI_MakeVertex(at(runs[n - 1], true)).Vertex();
    }
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for (size_t j = 0; j < n; ++j) {
        const bool back = runs[j].from > runs[j].to;
        BRepBuilderAPI_MakeEdge edge(runs[j].curve, back ? ends[j] : starts[j], back ? starts[j] : ends[j],
                                     std::min(runs[j].from, runs[j].to), std::max(runs[j].from, runs[j].to));
        if (!edge.IsDone())
            fail(op, no_solution, "a piece of the offset could not be built");
        builder.Add(wire, edge.Edge().Oriented(back ? TopAbs_REVERSED : TopAbs_FORWARD));
    }
    wire.Closed(closed);
    return wire;
}

// Whether an offset crosses itself, and how near it comes to its source, are measured on
// polylines: each edge walked so finely that the polyline is within `deflection` of it, in the
// plane's own two coordinates. OpenCascade's search for the least distance between an offset
// curve and a spline takes tens of milliseconds a pair.
struct Poly {
    std::vector<gp_Pnt2d> points;
    double low[2], high[2];
    int after = -1; // the polyline that follows this one in its wire, if any
};

Poly polyOf(const char *op, const TopoDS_Edge &edge, const gp_Pln &plane, double deflection) {
    BRepAdaptor_Curve curve(edge);
    // Steps that turn by no more than a tenth of a radian, beside the deflection: a walk by
    // deflection alone takes an S between two points for its chord.
    GCPnts_TangentialDeflection walk(curve, 0.1, deflection);
    if (walk.NbPoints() < 2)
        fail(op, no_solution, "a piece of the offset could not be walked along");
    Poly poly;
    for (int i = 1; i <= walk.NbPoints(); ++i) {
        double u, v;
        ElSLib::Parameters(plane, walk.Value(i), u, v);
        poly.points.emplace_back(u, v);
    }
    // In the direction the edge runs in its wire.
    if (edge.Orientation() == TopAbs_REVERSED)
        std::reverse(poly.points.begin(), poly.points.end());
    poly.low[0] = poly.low[1] = std::numeric_limits<double>::infinity();
    poly.high[0] = poly.high[1] = -std::numeric_limits<double>::infinity();
    for (const auto &p : poly.points) {
        poly.low[0] = std::min(poly.low[0], p.X());
        poly.low[1] = std::min(poly.low[1], p.Y());
        poly.high[0] = std::max(poly.high[0], p.X());
        poly.high[1] = std::max(poly.high[1], p.Y());
    }
    return poly;
}

// How far apart the boxes of two polylines are (0 where they overlap).
double boxGap(const Poly &a, const Poly &b) {
    const double dx = std::max({0.0, a.low[0] - b.high[0], b.low[0] - a.high[0]});
    const double dy = std::max({0.0, a.low[1] - b.high[1], b.low[1] - a.high[1]});
    return std::hypot(dx, dy);
}

double turn(const gp_Pnt2d &a, const gp_Pnt2d &b, const gp_Pnt2d &c) {
    return (b.X() - a.X()) * (c.Y() - a.Y()) - (b.Y() - a.Y()) * (c.X() - a.X());
}

// Whether two segments cross inside both: each has the other's ends strictly on its two sides.
bool segmentsCross(const gp_Pnt2d &a, const gp_Pnt2d &b, const gp_Pnt2d &c, const gp_Pnt2d &d) {
    const double d1 = turn(a, b, c), d2 = turn(a, b, d), d3 = turn(c, d, a), d4 = turn(c, d, b);
    // An end on the other segment's line, to the rounding of the products, is a touch.
    const double level = 1e-12 * a.Distance(b) * c.Distance(d);
    return d1 * d2 < 0 && d3 * d4 < 0 && std::min({std::abs(d1), std::abs(d2), std::abs(d3), std::abs(d4)}) > level;
}

double pointToSegment(const gp_Pnt2d &p, const gp_Pnt2d &a, const gp_Pnt2d &b) {
    const gp_Vec2d ab(a, b), ap(a, p);
    const double l2 = ab.SquareMagnitude();
    const double t = l2 > 0 ? std::min(1.0, std::max(0.0, ap.Dot(ab) / l2)) : 0;
    return p.Distance(a.Translated(ab * t));
}

// Whether two polylines cross. Of two that follow each other in a wire the two segments at
// their joint are passed over: their ends may lie 1e-7 apart there, which nearly parallel
// segments turn into a crossing.
bool cross(const Poly &a, const Poly &b, bool joined, bool self) {
    if (!self && boxGap(a, b) > 0)
        return false;
    const size_t na = a.points.size() - 1, nb = b.points.size() - 1;
    for (size_t i = 0; i < na; ++i)
        for (size_t k = self ? i + 2 : 0; k < nb; ++k) {
            if (joined && i + 1 == na && k == 0)
                continue;
            if (segmentsCross(a.points[i], a.points[i + 1], b.points[k], b.points[k + 1]))
                return true;
        }
    return false;
}

bool anyCross(const std::vector<Poly> &polys) {
    for (size_t i = 0; i < polys.size(); ++i)
        for (size_t k = i; k < polys.size(); ++k) {
            if (i == k ? cross(polys[i], polys[i], false, true)
                       : polys[i].after == int(k) ? cross(polys[i], polys[k], true, false)
                                                  : polys[k].after == int(i) ? cross(polys[k], polys[i], true, false)
                                                                             : cross(polys[i], polys[k], false, false))
                return true;
        }
    return false;
}

// The least distance between two sets of polylines.
double leastGap(const std::vector<Poly> &some, const std::vector<Poly> &others) {
    double least = std::numeric_limits<double>::infinity();
    for (const auto &a : some)
        for (const auto &b : others) {
            if (boxGap(a, b) >= least)
                continue;
            for (size_t i = 0; i + 1 < a.points.size(); ++i)
                for (size_t k = 0; k + 1 < b.points.size(); ++k) {
                    const auto &p = a.points[i], &q = a.points[i + 1], &r = b.points[k], &t = b.points[k + 1];
                    if (segmentsCross(p, q, r, t))
                        return 0;
                    least = std::min({least, pointToSegment(p, r, t), pointToSegment(q, r, t),
                                      pointToSegment(r, p, q), pointToSegment(t, p, q)});
                }
        }
    return least;
}

// The polylines of a wire's edges in the order it runs, each knowing the one after it.
void addPolys(const char *op, const TopoDS_Wire &wire, const gp_Pln &plane, double deflection,
              std::vector<Poly> &polys) {
    const size_t first = polys.size();
    for (BRepTools_WireExplorer it(wire); it.More(); it.Next())
        polys.push_back(polyOf(op, it.Current(), plane, deflection));
    for (size_t i = first; i + 1 < polys.size(); ++i)
        polys[i].after = int(i + 1);
    if (BRep_Tool::IsClosed(wire) && polys.size() > first + 1)
        polys.back().after = int(first);
}

// Where an end of a run goes when it is carried on, in the plane: along the line or circle it
// is of (`own`), along its tangent for any other curve.
struct Onward {
    Handle(Geom2d_Curve) curve;
    GeomAbs_CurveType type; // GeomAbs_Line, GeomAbs_Circle, or GeomAbs_OtherCurve for a tangent
    gp_Pnt from;            // the end
    gp_Vec out;             // the direction it leaves the end in
    double room;            // how far it may go
    double sign;            // of its curve's parameter in that direction
};

Onward onward(const Run &run, bool end, const gp_Pln &plane, double reach) {
    GeomAdaptor_Curve adaptor(run.curve);
    const gp_Pnt p = at(run, end);
    // Travel leaves the end forwards and the start backwards.
    const gp_Vec out = heading(run, end).Normalized() * (end ? 1.0 : -1.0);
    const double sign = (run.from > run.to ? -1.0 : 1.0) * (end ? 1.0 : -1.0);
    if (adaptor.GetType() == GeomAbs_Line)
        return {GeomAPI::To2d(run.curve, plane), GeomAbs_Line, p, out, reach, sign};
    if (adaptor.GetType() == GeomAbs_Circle)
        return {GeomAPI::To2d(run.curve, plane), GeomAbs_Circle, p, out,
                (2 * M_PI - std::abs(run.to - run.from)) * adaptor.Circle().Radius(), sign};
    Handle(Geom_Curve) line = new Geom_Line(p, gp_Dir(out));
    return {GeomAPI::To2d(line, plane), GeomAbs_OtherCurve, p, out, reach, 1.0};
}

// How far along `on` the point `meet` of it lies from the end, and the parameter of the run's
// own curve there; nothing if it lies behind the end or beyond the room. The parameters are
// read on the curve in space: the curve drawn in the plane need not keep them (a circle's
// start and sense follow the plane's axes).
std::optional<std::pair<double, double>> travel(const Run &run, bool end, const Onward &on, const gp_Pnt &meet) {
    const double u = end ? run.to : run.from;
    if (on.type == GeomAbs_Circle) {
        const gp_Circ circle = GeomAdaptor_Curve(run.curve).Circle();
        double turned = on.sign * (ElCLib::Parameter(circle, meet) - u);
        turned -= 2 * M_PI * std::floor(turned / (2 * M_PI));
        // A whole turn less the rounding is no turn.
        if (2 * M_PI - turned < tangent)
            turned = 0;
        const double length = turned * circle.Radius();
        if (length > on.room + same)
            return std::nullopt;
        return std::make_pair(length, u + on.sign * turned);
    }
    const double length = gp_Vec(on.from, meet).Dot(on.out);
    if (length < -same || length > on.room)
        return std::nullopt;
    return std::make_pair(std::max(0.0, length), u + on.sign * length);
}

// What an offset is measured against: the plane, how finely curves are walked, the source's
// polylines, and the least distance from the source that counts as the offset's distance.
struct Measure {
    gp_Pln plane;
    double deflection;
    std::vector<Poly> source;
    double least;
};

// The runs of one wire of the offset with their arcs replaced. `all` are the polylines of the
// whole offset, those of these runs from `first` on, and `added` the pieces put in so far: a
// piece that a run is carried on by must cross none of them, or the arc stays.
std::vector<Run> filled(const char *op, std::vector<Run> runs, bool closed, size_t first,
                        const std::vector<Poly> &all, std::vector<Poly> &added, const Measure &measure,
                        double distance, int gap) {
    const size_t n = runs.size();
    const gp_Pln &plane = measure.plane;
    const double reach = 2 * workSize + 100 * std::abs(distance);
    std::vector<std::vector<Run>> before(n), after(n); // lines put in front of and behind a run
    std::vector<bool> drop(n, false);
    for (size_t j = 0; j < n; ++j) {
        if (!runs[j].join)
            continue;
        if (gap == 1) {
            // The chord of half a circle or more would pass through the corner or behind it:
            // there the arc stays.
            if (std::abs(runs[j].to - runs[j].from) >= M_PI - tangent)
                continue;
            const gp_Pnt a = at(runs[j], false), b = at(runs[j], true);
            if (a.Distance(b) <= same) {
                drop[j] = true;
                continue;
            }
            Handle(Geom_Curve) chord = new Geom_Line(a, gp_Dir(gp_Vec(a, b)));
            runs[j] = {chord, 0, a.Distance(b), true};
            continue;
        }
        // Natural: where the neighbours cannot be carried on to a crossing, the arc stays, as
        // OpenCascade leaves it between lines and circles that do not cross.
        if (n < 3 || (!closed && (j == 0 || j + 1 == n)))
            continue;
        const size_t i = (j + n - 1) % n, k = (j + 1) % n;
        if (runs[i].join || runs[k].join)
            continue;
        // A neighbour that does not leave the arc in its direction was cut there by another
        // part of the offset: carried on, it would cross that part.
        if (heading(runs[i], true).Angle(heading(runs[j], false)) > tangent ||
            heading(runs[j], true).Angle(heading(runs[k], false)) > tangent)
            continue;
        const Onward a = onward(runs[i], true, plane, reach), b = onward(runs[k], false, plane, reach);
        Geom2dAPI_InterCurveCurve crossing(a.curve, b.curve, same);
        double best = std::numeric_limits<double>::infinity(), ua = 0, ub = 0, la = 0, lb = 0;
        gp_Pnt meet;
        for (int s = 1; s <= crossing.NbPoints(); ++s) {
            const gp_Pnt2d point = crossing.Point(s);
            const gp_Pnt here = ElSLib::Value(point.X(), point.Y(), plane);
            const auto onA = travel(runs[i], true, a, here), onB = travel(runs[k], false, b, here);
            if (!onA || !onB || onA->first + onB->first >= best)
                continue;
            best = onA->first + onB->first;
            la = onA->first;
            lb = onB->first;
            ua = onA->second;
            ub = onB->second;
            meet = here;
        }
        if (!std::isfinite(best))
            continue;
        // What an end is carried on by: nothing where it is at the crossing already.
        std::optional<Run> onA, onB;
        if (la > same) {
            Handle(Geom_Curve) line = new Geom_Line(a.from, gp_Dir(a.out));
            if (a.type == GeomAbs_OtherCurve)
                onA = Run{line, 0, la, false};
            else
                onA = Run{runs[i].curve, runs[i].to, ua, false};
        }
        if (lb > same) {
            Handle(Geom_Curve) line = new Geom_Line(meet, gp_Dir(-b.out));
            if (b.type == GeomAbs_OtherCurve)
                onB = Run{line, 0, lb, false};
            else
                onB = Run{runs[k].curve, ub, runs[k].from, false};
        }
        // Carried on, a piece must stay at the distance from the source or beyond it, and
        // cross nothing else of the offset.
        std::vector<Poly> pieces;
        bool sound = true;
        for (const auto *piece : {&onA, &onB}) {
            if (!*piece)
                continue;
            pieces.push_back(polyOf(op, edgeOf(op, **piece), plane, measure.deflection));
            if (leastGap({pieces.back()}, measure.source) < measure.least)
                sound = false;
            for (size_t other = 0; other < all.size() && sound; ++other)
                if (other != first + i && other != first + j && other != first + k &&
                    cross(pieces.back(), all[other], false, false))
                    sound = false;
            for (const auto &other : added)
                if (sound && cross(pieces.back(), other, false, false))
                    sound = false;
        }
        if (!sound)
            continue;
        added.insert(added.end(), pieces.begin(), pieces.end());
        if (onA && a.type != GeomAbs_OtherCurve)
            runs[i].to = ua;
        else if (onA)
            after[i].push_back(*onA);
        if (onB && b.type != GeomAbs_OtherCurve)
            runs[k].from = ub;
        else if (onB)
            before[k].push_back(*onB);
        drop[j] = true;
    }
    std::vector<Run> result;
    for (size_t j = 0; j < n; ++j) {
        if (drop[j])
            continue;
        result.insert(result.end(), before[j].begin(), before[j].end());
        result.push_back(runs[j]);
        result.insert(result.end(), after[j].begin(), after[j].end());
    }
    return result;
}

bool linesAndCircles(const TopoDS_Shape &shape) {
    for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next()) {
        const auto type = BRepAdaptor_Curve(TopoDS::Edge(it.Current())).GetType();
        if (type != GeomAbs_Line && type != GeomAbs_Circle)
            return false;
    }
    return true;
}

// Whether some point inside a closed curve (its polylines) lies deeper in it than `depth`, by
// more than the share that counts: whether anything is left of the curve offset inwards by
// that much. A search, not a sample: the depth of a point changes no faster than the point
// moves, so a square whose middle is too shallow by more than its half diagonal holds no such
// point and is left; any other is quartered, down to squares smaller than what counts.
bool deeperThan(const std::vector<Poly> &loop, double depth) {
    const double wanted = depth * (1 + nearShare);
    double low[2] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    double high[2] = {-low[0], -low[1]};
    for (const auto &poly : loop)
        for (int k = 0; k < 2; ++k) {
            low[k] = std::min(low[k], poly.low[k]);
            high[k] = std::max(high[k], poly.high[k]);
        }
    // How deep inside a point is; negative outside.
    auto depthAt = [&](const gp_Pnt2d &p) {
        double nearest = std::numeric_limits<double>::infinity();
        bool inside = false;
        for (const auto &poly : loop)
            for (size_t i = 0; i + 1 < poly.points.size(); ++i) {
                const auto &a = poly.points[i], &b = poly.points[i + 1];
                nearest = std::min(nearest, pointToSegment(p, a, b));
                // A ray towards +x crosses the segment.
                if ((a.Y() > p.Y()) != (b.Y() > p.Y()) &&
                    p.X() < a.X() + (b.X() - a.X()) * (p.Y() - a.Y()) / (b.Y() - a.Y()))
                    inside = !inside;
            }
        return inside ? nearest : -nearest;
    };
    struct Square {
        double x, y, half;
    };
    const double half = std::max(high[0] - low[0], high[1] - low[1]) / 2;
    std::vector<Square> open{{(low[0] + high[0]) / 2, (low[1] + high[1]) / 2, half}};
    while (!open.empty()) {
        const Square square = open.back();
        open.pop_back();
        const double here = depthAt(gp_Pnt2d(square.x, square.y));
        if (here > wanted)
            return true;
        if (here + square.half * M_SQRT2 <= wanted || square.half < depth * nearShare / 2)
            continue;
        for (const double dx : {-0.5, 0.5})
            for (const double dy : {-0.5, 0.5})
                open.push_back({square.x + dx * square.half, square.y + dy * square.half, square.half / 2});
    }
    return false;
}

// Twice the area a closed chain of runs goes round, from a polygon through 32 points of each
// (exact for straight runs: this tells a loop of two lines lying on each other from a real one).
double twiceArea(const std::vector<Run> &runs, const gp_Pln &plane) {
    std::vector<gp_Pnt2d> points;
    for (const auto &run : runs)
        for (int s = 0; s < 32; ++s) {
            double u, v;
            ElSLib::Parameters(plane, run.curve->Value(run.from + (run.to - run.from) * s / 32.0), u, v);
            points.emplace_back(u, v);
        }
    double sum = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto &p = points[i], &q = points[(i + 1) % points.size()];
        sum += p.X() * q.Y() - q.X() * p.Y();
    }
    return sum;
}

// The offset of `source` (a wire, or a planar face for its boundary; `outline` is the curve
// itself) with its gaps closed the way `gap` says, as wires: none where nothing is left of a
// closed curve. What is returned crosses itself nowhere and is nowhere nearer to the curve than
// the distance: OpenCascade's offset is checked for both, since its own check passes either.
// `enlarged`: the size OpenCascade is given the curve at, as a multiple of the working size.
std::vector<TopoDS_Wire> offsetWiresAt(const char *op, const TopoDS_Shape &source, const TopoDS_Shape &outline,
                                       const gp_Pln &plane, bool open, double distance, int gap, double enlarged) {
    const bool own = gap == 1 || (gap == 2 && !linesAndCircles(source));
    const auto join = gap == 2 && !own ? GeomAbs_Intersection : GeomAbs_Arc;
    gp_Trsf enlarge;
    enlarge.SetScale(gp_Pnt(0, 0, 0), enlarged);
    const TopoDS_Shape given = enlarged == 1 ? source : BRepBuilderAPI_Transform(source, enlarge, true).Shape();
    BRepOffsetAPI_MakeOffset offset;
    if (given.ShapeType() == TopAbs_FACE)
        offset.Init(TopoDS::Face(given), join, open);
    else {
        offset.Init(join, open);
        offset.AddWire(TopoDS::Wire(given));
    }
    offset.Perform(distance * enlarged);
    if (!offset.IsDone() || offset.Shape().IsNull())
        fail(op, no_solution, "the curve could not be offset");
    std::optional<BRepBuilderAPI_Transform> reduce;
    if (enlarged != 1)
        reduce.emplace(offset.Shape(), enlarge.Inverted(), true);
    const TopoDS_Shape made = reduce ? reduce->Shape() : offset.Shape();
    // The arcs that close the corners: what the vertices of the source generated, each one
    // checked to be a circle of the offset's radius about its vertex.
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> joins;
    if (own)
        for (TopExp_Explorer v(given, TopAbs_VERTEX); v.More(); v.Next()) {
            const gp_Pnt corner = BRep_Tool::Pnt(TopoDS::Vertex(v.Current()));
            for (const auto &piece : offset.Generated(v.Current())) {
                if (piece.ShapeType() != TopAbs_EDGE)
                    continue;
                BRepAdaptor_Curve arc(TopoDS::Edge(piece));
                if (arc.GetType() != GeomAbs_Circle ||
                    std::abs(arc.Circle().Radius() - std::abs(distance) * enlarged) > same * enlarged ||
                    arc.Circle().Location().Distance(corner) > same * enlarged)
                    fail(op, no_solution, "a corner of the offset is not closed by an arc about it");
                joins.Add(reduce ? reduce->ModifiedShape(piece) : piece);
            }
        }
    // Walked to a tenth of what "nearer than the distance" means, within what a curve of the
    // working size can be walked to in a few thousand steps.
    Measure measure{plane, std::min(1e-3, std::max(1e-6, 0.1 * nearShare * std::abs(distance))), {}, 0};
    measure.least = std::abs(distance) * (1 - nearShare) - 2 * measure.deflection;
    for (TopExp_Explorer w(outline, TopAbs_WIRE); w.More(); w.Next())
        addPolys(op, TopoDS::Wire(w.Current()), plane, measure.deflection, measure.source);
    struct Loop {
        TopoDS_Wire wire;
        std::vector<Run> runs;
        size_t first;
        bool closed;
    };
    std::vector<Loop> loops;
    std::vector<Poly> all;
    for (TopExp_Explorer w(made, TopAbs_WIRE); w.More(); w.Next()) {
        Loop loop{TopoDS::Wire(w.Current()), runsOf(op, TopoDS::Wire(w.Current()), joins), all.size(), false};
        loop.closed = BRep_Tool::IsClosed(loop.wire);
        // Where two parts of the offset just touch, OpenCascade leaves a loop of two lines that
        // lie on each other: it goes round nothing.
        if (!open && loop.closed && std::abs(twiceArea(loop.runs, plane)) <= 1e-9 * workSize * workSize)
            continue;
        addPolys(op, loop.wire, plane, measure.deflection, all);
        loops.push_back(loop);
    }
    if (loops.empty())
        return {};
    if (anyCross(all))
        fail(op, no_solution, "the offset crosses itself at this distance");
    if (leastGap(all, measure.source) < measure.least)
        fail(op, no_solution, "the offset comes nearer to the curve than the distance: it would cross itself");
    std::vector<TopoDS_Wire> wires;
    if (!own) {
        for (const auto &loop : loops)
            wires.push_back(loop.wire);
        return wires;
    }
    std::vector<Poly> added, final;
    for (const auto &loop : loops) {
        const auto runs = filled(op, loop.runs, loop.closed, loop.first, all, added, measure, distance, gap);
        if (runs.empty())
            continue;
        wires.push_back(wireOf(op, runs, loop.closed, "the pieces of the offset do not join up"));
        if (!BRepCheck_Analyzer(wires.back()).IsValid())
            fail(op, no_solution, "the offset with its gaps closed is not a valid curve");
        addPolys(op, wires.back(), plane, measure.deflection, final);
    }
    if (anyCross(final))
        fail(op, no_solution, "the offset with its gaps closed crosses itself at this distance");
    return wires;
}

// OpenCascade's offset gives up on, or returns crossing itself, some curves at some sizes and
// not at others (an open line, arc and spline at the working size; a spline blob outwards at
// that size and not at three times it). It is asked at three sizes; the first answer that
// passes the checks is taken. Where none does, nothing is left of a closed curve offset
// inwards if no point inside it is deeper than the distance: OpenCascade's returning nothing,
// or giving up, does not tell (it returned nothing for a spline blob with most of it left).
std::vector<TopoDS_Wire> offsetWires(const char *op, const TopoDS_Shape &source, const TopoDS_Shape &outline,
                                     const gp_Pln &plane, bool open, double distance, int gap) {
    std::optional<Failure> refused;
    for (const double enlarged : {1.0, 3.0, 1 / 3.0}) {
        try {
            auto wires = offsetWiresAt(op, source, outline, plane, open, distance, gap, enlarged);
            if (!wires.empty())
                return wires;
        } catch (const Failure &failure) {
            if (!refused)
                refused = failure;
        }
    }
    if (!open && distance < 0) {
        std::vector<Poly> loop;
        for (TopExp_Explorer w(outline, TopAbs_WIRE); w.More(); w.Next())
            addPolys(op, TopoDS::Wire(w.Current()), plane, 1e-5 * workSize, loop);
        if (!deeperThan(loop, -distance))
            return {};
    }
    if (refused)
        throw *refused;
    fail(op, no_solution, "the curve could not be offset");
}

// The copy of a shape that is worked on: moved so that its middle is the origin and scaled to
// `workSize`. `factor` is the scale, `back` takes a result to where the shape is.
TopoDS_Shape atWorkSize(const char *op, const TopoDS_Shape &shape, double &factor, gp_Trsf &back) {
    // The box of the curves themselves, without the shape's tolerances.
    Bnd_Box box;
    BRepBndLib::AddOptimal(shape, box, false, false);
    if (box.IsVoid() || box.SquareExtent() <= 0)
        fail(op, invalid_input, "the curve has no size");
    double x0, y0, z0, x1, y1, z1;
    box.Get(x0, y0, z0, x1, y1, z1);
    const gp_Pnt middle((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
    factor = workSize / std::sqrt(box.SquareExtent());
    gp_Trsf scale, move;
    scale.SetScale(middle, factor);
    move.SetTranslation(middle, gp_Pnt(0, 0, 0));
    const gp_Trsf there = move * scale;
    back = there.Inverted();
    return BRepBuilderAPI_Transform(shape, there, true).Shape();
}

// A wire of the copy drawn afresh from its curves. Scaling leaves two things as they were,
// and OpenCascade's offset answers differently by both: the tolerances of vertices and edges
// (scaled with the copy), and the knots of a spline (not scaled: a spline through points a
// thousand times as far apart has knots a thousand times as far apart, at any size it is then
// scaled to). The knots are scaled like the copy, and the edges and vertices made new.
TopoDS_Wire redrawn(const char *op, const TopoDS_Wire &wire, double factor) {
    auto runs = runsOf(op, wire, {});
    for (auto &run : runs) {
        Handle(Geom_Curve) curve = run.curve;
        while (auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve))
            curve = trimmed->BasisCurve();
        const auto spline = Handle(Geom_BSplineCurve)::DownCast(curve);
        if (spline.IsNull())
            continue;
        const auto scaled = Handle(Geom_BSplineCurve)::DownCast(spline->Copy());
        NCollection_Array1<double> knots = scaled->Knots();
        for (int i = knots.Lower(); i <= knots.Upper(); ++i)
            knots(i) *= factor;
        scaled->SetKnots(knots);
        run = {Handle(Geom_Curve)(scaled), run.from * factor, run.to * factor, false};
    }
    return wireOf(op, runs, BRep_Tool::IsClosed(wire), "the pieces of the curve do not join up");
}

// Which side of an open curve its offset lies on, seen against the plane's own normal: 1 to
// the right of the way the curve runs, -1 to the left, 0 if no edge of the offset has its
// middle nearest to the inside of an edge of the curve.
int sideOf(const char *op, const std::vector<TopoDS_Wire> &wires, const TopoDS_Wire &curve, const gp_Pln &plane) {
    const double deflection = 1e-5 * workSize;
    std::vector<Poly> source, made;
    addPolys(op, curve, plane, deflection, source);
    for (const auto &wire : wires)
        addPolys(op, wire, plane, deflection, made);
    for (const auto &poly : made) {
        const gp_Pnt2d q = poly.points.size() > 2 ? poly.points[poly.points.size() / 2]
                                                  : gp_Pnt2d((poly.points[0].XY() + poly.points[1].XY()) / 2);
        double least = std::numeric_limits<double>::infinity(), aside = 0;
        for (const auto &from : source)
            for (size_t i = 0; i + 1 < from.points.size(); ++i) {
                const auto &a = from.points[i], &b = from.points[i + 1];
                const double gap = pointToSegment(q, a, b);
                if (gap >= least)
                    continue;
                least = gap;
                const gp_Vec2d ab(a, b), aq(a, q);
                const double t = aq.Dot(ab) / ab.SquareMagnitude();
                // Beside an end of an edge of the curve the nearest point is a corner or an
                // end, which has no one direction.
                const bool inside = (i > 0 || t > 0) && (i + 2 < from.points.size() || t < 1);
                aside = inside ? turn(a, b, q) / ab.Magnitude() : 0;
            }
        if (std::abs(aside) > same)
            return aside < 0 ? 1 : -1;
    }
    return 0;
}
} // namespace

// A planar curve at a distance from itself. One that is open keeps its ends: OpenCascade
// otherwise offsets it to both sides at once and closes the two round the ends, whatever the
// sign. For a closed curve positive is outwards. For an open one positive is to the right of
// the way it runs, seen against `up` (from where `up` points to); where `up` is zero or lies in
// the curve's plane, against that normal of its plane whose largest component is positive.
// OpenCascade picks the side by a plane it finds for the curve, whose normal turns over as the
// curve's plane is tilted in the world: the side is read off the result and the offset made
// again with the other sign where it is the wrong one.
static EvolutionData offsetCurveSeen(const char *op, OcctKernel &kernel, uint32_t wireId, double distance,
                                     int gap, Vec3 up) {
    const auto &shape = get(kernel, wireId);
    if (shape.ShapeType() != TopAbs_WIRE || !std::isfinite(distance) || distance == 0 || gap < 0 || gap > 2 ||
        !finite(up))
        fail(op, invalid_input, "a planar curve and a nonzero distance are required");
    double factor;
    gp_Trsf back;
    // The copy also keeps the curve given as it is: OpenCascade marks a wire it has put into a
    // face as no longer free.
    const auto wire = redrawn(op, TopoDS::Wire(atWorkSize(op, shape, factor, back)), factor);
    const bool open = !BRep_Tool::IsClosed(wire);
    BRepBuilderAPI_FindPlane found(wire);
    if (!found.Found())
        fail(op, unsupported_geometry, "the curve is not planar");
    const gp_Pln plane = found.Plane()->Pln();
    auto wires = offsetWires(op, wire, wire, plane, open, distance * factor, gap);
    if (open && !wires.empty()) {
        gp_Vec normal(plane.Axis().Direction());
        const gp_Vec given(up.x, up.y, up.z);
        const double along = given.Magnitude() > 0 ? normal.Dot(given) / given.Magnitude() : 0;
        if (std::abs(along) > tangent) {
            if (along < 0)
                normal.Reverse();
        } else {
            const double largest = std::max({std::abs(normal.X()), std::abs(normal.Y()), std::abs(normal.Z())});
            const double first = std::abs(normal.X()) == largest   ? normal.X()
                                 : std::abs(normal.Y()) == largest ? normal.Y()
                                                                   : normal.Z();
            if (first < 0)
                normal.Reverse();
        }
        // Seen against the plane's own normal, and the other way round if that is not the one.
        const int side = sideOf(op, wires, wire, plane) * (normal.Dot(gp_Vec(plane.Axis().Direction())) > 0 ? 1 : -1);
        if (side * distance < 0)
            wires = offsetWires(op, wire, wire, plane, open, -distance * factor, gap);
    }
    if (wires.empty())
        fail(op, no_solution, "nothing is left of the curve at this distance");
    BRep_Builder builder;
    TopoDS_Compound all;
    builder.MakeCompound(all);
    for (const auto &made : wires)
        builder.Add(all, made);
    return stored(kernel, BRepBuilderAPI_Transform(all, back, true).Shape());
}

EvolutionData offsetCurveOn(OcctKernel &kernel, uint32_t wireId, double distance, int gap, Vec3 up) {
    return offsetCurveSeen("offsetCurveOn", kernel, wireId, distance, gap, up);
}

EvolutionData offsetCurve(OcctKernel &kernel, uint32_t wireId, double distance, int gap) {
    return offsetCurveSeen("offsetCurve", kernel, wireId, distance, gap, Vec3{0, 0, 0});
}

// A planar region grown (positive) or shrunk by a distance: its outer boundary offset
// outwards and each hole inwards, or the other way. A boundary may vanish (a hole that closes,
// a region of which nothing is left) or come apart (a region pinched in two), and boundaries
// may run into each other (a hole that reaches the outer boundary or another hole): then the
// result is what is inside the outer loops and outside every hole's loop, by a Boolean of the
// faces. OpenCascade's check of a face does not see wires that cross, so whether they do is
// measured between every two loops.
EvolutionData offsetRegion(OcctKernel &kernel, uint32_t faceId, double distance, int gap, double tolerance) {
    const char *op = "offsetRegion";
    const auto &shape = get(kernel, faceId);
    if (shape.ShapeType() != TopAbs_FACE || !std::isfinite(distance) || std::abs(distance) <= tolerance ||
        gap < 0 || gap > 2)
        fail(op, invalid_input, "a planar region and nonzero distance are required");
    if (BRepAdaptor_Surface(TopoDS::Face(shape)).GetType() != GeomAbs_Plane)
        fail(op, invalid_input, "region must be planar");
    double factor;
    gp_Trsf back;
    const auto input = TopoDS::Face(atWorkSize(op, shape, factor, back).Oriented(TopAbs_FORWARD));
    const gp_Pln plane = BRepAdaptor_Surface(input).Plane();
    const double by = distance * factor;
    const auto outer = BRepTools::OuterWire(input);
    // The loops a boundary becomes, each as the face it goes round.
    auto loopsOf = [&](const TopoDS_Wire &boundary, bool hole) {
        // Every boundary is offset as a disk; its role alone chooses the sign.
        const TopoDS_Wire wire = redrawn(op, boundary, factor);
        BRepBuilderAPI_MakeFace disk(plane, wire, true);
        std::vector<TopoDS_Face> faces;
        for (const auto &made : offsetWires(op, disk.Face(), wire, plane, false, hole ? -by : by, gap))
            faces.push_back(BRepBuilderAPI_MakeFace(plane, made, true).Face());
        return faces;
    };
    const auto outers = loopsOf(outer, false);
    if (outers.empty())
        fail(op, no_solution, "nothing is left of the region at this distance");
    std::vector<TopoDS_Face> holes;
    for (TopExp_Explorer it(input, TopAbs_WIRE); it.More(); it.Next())
        if (!it.Current().IsSame(outer))
            for (const auto &hole : loopsOf(TopoDS::Wire(it.Current()), true))
                holes.push_back(hole);
    // Loops that keep clear of each other bound one face as they are, if each hole lies inside
    // the outer loop and in no other hole: one point of a loop tells, since they do not meet.
    auto inside = [&](const TopoDS_Face &loop, const TopoDS_Face &of) {
        TopExp_Explorer vertex(loop, TopAbs_VERTEX);
        return BRepClass_FaceClassifier(of, BRep_Tool::Pnt(TopoDS::Vertex(vertex.Current())), same).State() ==
               TopAbs_IN;
    };
    const double deflection = 1e-5 * workSize;
    auto walked = [&](const TopoDS_Face &loop) {
        std::vector<Poly> polys;
        addPolys(op, BRepTools::OuterWire(loop), plane, deflection, polys);
        return polys;
    };
    bool clear = outers.size() == 1;
    std::vector<std::vector<Poly>> holePolys;
    for (const auto &hole : holes)
        holePolys.push_back(walked(hole));
    const auto outerPolys = walked(outers[0]);
    // Polylines within the deflection of two loops are this far apart where the loops meet.
    const double apart = 4 * deflection;
    for (size_t i = 0; i < holes.size() && clear; ++i) {
        clear = leastGap(holePolys[i], outerPolys) > apart && inside(holes[i], outers[0]);
        for (size_t k = 0; k < holes.size() && clear; ++k)
            if (k != i)
                clear = leastGap(holePolys[i], holePolys[k]) > apart && !inside(holes[i], holes[k]);
    }
    GProp_GProps before, after;
    TopoDS_Shape result;
    BRep_Builder builder;
    if (clear) {
        BRepBuilderAPI_MakeFace face(plane, BRepTools::OuterWire(outers[0]), true);
        for (const auto &hole : holes)
            face.Add(TopoDS::Wire(BRepTools::OuterWire(hole).Reversed()));
        result = face.Face();
    } else {
        TopoDS_Compound all;
        builder.MakeCompound(all);
        for (const auto &face : outers)
            builder.Add(all, face);
        result = all;
        for (const auto &hole : holes) {
            BRepAlgoAPI_Cut cut(result, hole);
            if (!cut.IsDone())
                fail(op, no_solution, "the region's boundaries run into each other and could not be cut");
            result = cut.Shape();
        }
        ShapeUpgrade_UnifySameDomain unify(result, true, true, false);
        unify.Build();
        result = unify.Shape();
        std::vector<TopoDS_Shape> faces;
        for (TopExp_Explorer it(result, TopAbs_FACE); it.More(); it.Next())
            faces.push_back(it.Current());
        if (faces.empty())
            fail(op, no_solution, "nothing is left of the region at this distance");
        if (faces.size() == 1)
            result = faces[0];
        else {
            TopoDS_Compound parts;
            builder.MakeCompound(parts);
            for (const auto &face : faces)
                builder.Add(parts, face);
            result = parts;
        }
        // What the cut leaves is no more than the outer loops go round and no less than that
        // without all that the holes' loops go round.
        double most = 0;
        for (const auto &face : outers) {
            GProp_GProps area;
            BRepGProp::SurfaceProperties(face, area);
            most += area.Mass();
        }
        double least = most;
        for (const auto &hole : holes) {
            GProp_GProps area;
            BRepGProp::SurfaceProperties(hole, area);
            least -= area.Mass();
        }
        BRepGProp::SurfaceProperties(result, after);
        const double slack = 1e-9 * workSize * workSize;
        if (after.Mass() > most + slack || after.Mass() < least - slack)
            fail(op, no_solution, "the region's boundaries run into each other and could not be cut");
    }
    BRepGProp::SurfaceProperties(input, before);
    after = GProp_GProps();
    BRepGProp::SurfaceProperties(result, after);
    if (after.Mass() <= 0 || (after.Mass() - before.Mass()) * distance <= 0)
        fail(op, no_solution, "the offset region does not grow or shrink as asked");
    result = BRepBuilderAPI_Transform(result, back, true).Shape();
    if (!BRepCheck_Analyzer(result).IsValid()) {
        // A face made of an offset's wires fails OpenCascade's check where an edge's curve in
        // the plane is not within the edge's tolerance of its curve in space; its repair
        // settles that.
        ShapeFix_Shape fix(result);
        fix.Perform();
        result = fix.Shape();
        if (!BRepCheck_Analyzer(result).IsValid())
            fail(op, no_solution, "offset region is invalid");
    }
    return stored(kernel, result);
}
EvolutionData offsetVertex(OcctKernel &kernel, uint32_t wireId, uint32_t vertexId, double distance,
                           double tolerance) {
    const char *op = "offsetVertex";
    const auto &shape = get(kernel, wireId);
    if (shape.ShapeType() != TopAbs_WIRE || !std::isfinite(distance) || distance <= tolerance)
        fail(op, invalid_input, "an open or closed wire and positive distance are required");
    const auto vertices = subShapes(kernel, op, shape, {vertexId}, TopAbs_VERTEX);
    const auto point = BRep_Tool::Pnt(TopoDS::Vertex(vertices[0]));
    // The edges in the wire's order, each with its length, and which of them start at the vertex.
    struct Run {
        TopoDS_Edge edge;
        double length;
        gp_Pnt from, to;
    };
    std::vector<Run> runs;
    for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next()) {
        const auto edge = it.Current();
        BRepAdaptor_Curve curve(edge);
        const double a = curve.FirstParameter(), b = curve.LastParameter();
        const bool reversed = edge.Orientation() == TopAbs_REVERSED;
        runs.push_back({edge, GCPnts_AbscissaPoint::Length(curve, a, b, tolerance), curve.Value(reversed ? b : a),
                        curve.Value(reversed ? a : b)});
    }
    const int count = int(runs.size());
    const bool closed = count > 1 && runs.front().from.Distance(runs.back().to) <= tolerance;
    int after = -1;
    for (int i = 0; i < count; ++i)
        if (runs[size_t(i)].from.Distance(point) <= tolerance) {
            if (runs[size_t(i)].to.Distance(point) <= tolerance)
                fail(op, invalid_input, "a closed single edge has no two-sided corner");
            after = i;
        }
    const bool atJoint = after > 0 || (after == 0 && closed);
    if (!atJoint)
        fail(op, invalid_input, "select a vertex between exactly two segments");
    // Walks on along the wire from the vertex, past as many vertices as the distance takes: which
    // edge it ends in, and how far along that edge in the direction of the walk.
    struct Stop {
        int run;
        double along;
    };
    const auto walk = [&](int step) {
        double left = distance;
        int i = step > 0 ? after : after - 1;
        for (int walked = 0; walked < count; ++walked) {
            const int at = closed ? (i + count) % count : i;
            if (at < 0 || at >= count)
                fail(op, invalid_input, "the distance runs past the end of the curve");
            const double length = runs[size_t(at)].length;
            if (left < length - tolerance)
                return Stop{at, left};
            left -= length;
            // Landing on a vertex needs no new one.
            if (left <= tolerance)
                return Stop{at, length};
            i += step;
        }
        fail(op, invalid_input, "the distance goes all the way round the curve");
    };
    const Stop onward = walk(1), back = walk(-1);
    double round = 0;
    for (const auto &run : runs)
        round += run.length;
    if (closed && 2 * distance >= round - tolerance)
        fail(op, invalid_input, "the two sides meet on the far side of the curve");
    // An edge cut in two, in the order of the wire, at a length from its start along the wire.
    const auto cutAt = [&](const Run &run, std::vector<double> lengths, BRepBuilderAPI_MakeWire &output) {
        BRepAdaptor_Curve curve(run.edge);
        const double a = curve.FirstParameter(), b = curve.LastParameter();
        const bool reversed = run.edge.Orientation() == TopAbs_REVERSED;
        std::sort(lengths.begin(), lengths.end());
        double lo, hi;
        const auto basis = BRep_Tool::Curve(run.edge, lo, hi);
        if (basis.IsNull())
            fail(op, invalid_input, "edge has no 3D curve");
        double from = reversed ? b : a;
        for (size_t k = 0; k <= lengths.size(); ++k) {
            double to = reversed ? a : b;
            if (k < lengths.size()) {
                GCPnts_AbscissaPoint locate(tolerance, curve, reversed ? -lengths[k] : lengths[k], reversed ? b : a);
                if (!locate.IsDone())
                    fail(op, no_solution, "arc-length split failed");
                to = locate.Parameter();
            }
            auto piece = BRepBuilderAPI_MakeEdge(basis, std::min(from, to), std::max(from, to)).Edge();
            output.Add(reversed ? TopoDS::Edge(piece.Reversed()) : piece);
            from = to;
        }
    };
    BRepBuilderAPI_MakeWire output;
    for (int i = 0; i < count; ++i) {
        const Run &run = runs[size_t(i)];
        std::vector<double> cuts;
        if (i == onward.run && onward.along < run.length - tolerance)
            cuts.push_back(onward.along);
        if (i == back.run && back.along < run.length - tolerance)
            cuts.push_back(run.length - back.along);
        if (cuts.empty())
            output.Add(run.edge);
        else
            cutAt(run, cuts, output);
    }
    if (!output.IsDone() || !BRepCheck_Analyzer(output.Wire()).IsValid())
        fail(op, no_solution, "split wire is invalid");
    return stored(kernel, output.Wire());
}
} // namespace BRepExt
