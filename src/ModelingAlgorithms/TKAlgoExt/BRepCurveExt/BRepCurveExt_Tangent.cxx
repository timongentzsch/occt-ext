// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_PlanarCurve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Tool.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dGcc_Circ2d2TanOn.hxx>
#include <Geom2dGcc_Circ2d2TanRad.hxx>
#include <Geom2dGcc_Circ2d3Tan.hxx>
#include <Geom2dGcc_QualifiedCurve.hxx>
#include <Geom2d_CartesianPoint.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <gp_Circ.hxx>
#include <gp_Circ2d.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <limits>
#include <vector>

namespace BRepExt {
namespace {
// Everything here is worked out in a frame of the curves' own (`Flat`): the plane as z = 0 and
// the middle of the curves at the origin, so that the answer does not depend on where the
// model sits. A point of a curve is known to about a double's resolution of the curves' size
// (`extent`), and no better, so the tolerances are shares of the radius or of that size, none a
// length of its own:
// - Newton stops where the centre is at the radius from each curve and square to it there to
//   1e-13 of the extent;
// - a circle is accepted where that holds to 1e-9 of its radius (what a caller may hold the
//   result to) plus 1e-12 of the extent.
constexpr double exactShare = 1e-9, acceptedShare = 1e-12, resolvedShare = 1e-13;

struct Flat {
    gp_Trsf toLocal;
    std::vector<Handle(Geom_Curve)> curves; // in the frame; what is across the plane is rounding
    std::vector<double> first, last;
    double extent = 0;
};

struct Circle {
    double x, y, radius; // in the frame
};

struct At {
    double x, y, dx, dy, ddx, ddy;
};
At at(const Flat &flat, int i, double t) {
    gp_Pnt p;
    gp_Vec d1, d2;
    flat.curves[i]->D2(t, p, d1, d2);
    return At{p.X(), p.Y(), d1.X(), d1.Y(), d2.X(), d2.Y()};
}

// Newton's method on the centre, the radius and the three parameters, for the centre to be at
// the radius from each curve and square to it there. Where two of the curves touch each other
// the system is singular at the circle through that point and the steps only halve: they are
// followed until they stop shrinking, not only until the equations hold, since there the
// equations hold long before the centre is in place.
bool settle(const Flat &flat, double t[3], Circle &circle) {
    double &cx = circle.x, &cy = circle.y, &r = circle.radius;
    double before = std::numeric_limits<double>::infinity();
    bool still = false;
    for (int step = 0; step < 80; ++step) {
        math_Matrix m(0, 5, 0, 5, 0.);
        math_Vector rhs(0, 5), x(0, 5);
        double residual = 0;
        for (int c = 0; c < 3; ++c) {
            const At p = at(flat, c, t[c]);
            const double rx = cx - p.x, ry = cy - p.y;
            const double square = rx * p.dx + ry * p.dy, away = rx * rx + ry * ry - r * r;
            const double speed = std::hypot(p.dx, p.dy);
            residual = std::max({residual, std::abs(square) / std::max(speed, 1e-300),
                                 std::abs(std::sqrt(rx * rx + ry * ry) - r)});
            m(2 * c, 0) = p.dx;
            m(2 * c, 1) = p.dy;
            m(2 * c, 3 + c) = -speed * speed + rx * p.ddx + ry * p.ddy;
            rhs(2 * c) = -square;
            m(2 * c + 1, 0) = 2 * rx;
            m(2 * c + 1, 1) = 2 * ry;
            m(2 * c + 1, 2) = -2 * r;
            m(2 * c + 1, 3 + c) = -2 * square;
            rhs(2 * c + 1) = -away;
        }
        const bool met = residual <= resolvedShare * flat.extent;
        if (met && still)
            return true;
        const math_Gauss solved(m, 1e-300);
        if (!solved.IsDone())
            return met;
        solved.Solve(rhs, x);
        const double moved = std::hypot(x(0), x(1)) + std::abs(x(2));
        cx += x(0);
        cy += x(1);
        r += x(2);
        for (int c = 0; c < 3; ++c) {
            t[c] += x(3 + c);
            if (!(t[c] >= flat.first[c] && t[c] <= flat.last[c]))
                return false;
        }
        if (!(r > 0) || !std::isfinite(r))
            return false;
        still = moved <= 1e-15 * flat.extent || moved >= before;
        before = moved;
    }
    return false;
}

// The circle through three points; false where they lie on a line.
bool through(const At p[3], Circle &circle) {
    const double ax = p[1].x - p[0].x, ay = p[1].y - p[0].y, bx = p[2].x - p[0].x, by = p[2].y - p[0].y;
    const double d = 2 * (ax * by - ay * bx);
    if (std::abs(d) <= 1e-12 * std::hypot(ax, ay) * std::hypot(bx, by))
        return false;
    circle.x = p[0].x + (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    circle.y = p[0].y + (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    circle.radius = std::hypot(circle.x - p[0].x, circle.y - p[0].y);
    return true;
}

// How many places of a curve to look at: four to each span of a spline (a span is one
// polynomial, which turns only so often), and no fewer than 64.
int placesOn(const Handle(Geom_Curve) & curve) {
    Handle(Geom_Curve) basis = curve;
    if (const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve))
        basis = trimmed->BasisCurve();
    int spans = 16;
    if (const auto spline = Handle(Geom_BSplineCurve)::DownCast(basis))
        spans = std::max(spans, spline->NbKnots() - 1);
    return std::min(4 * spans, 2048);
}

// The circles tangent to three curves of a plane. OpenCascade's Geom2dGcc_Circ2d3Tan iterates
// from one start on each curve and finds no circle from most starts once a spline is among
// them. Two searches of our own give the starts for `settle`:
//
// 1. A grid of 8 places on each curve, and the circle through each three of them. This finds
//    circles of every kind, those that cut through a curve elsewhere among them.
//
// 2. The circles that cut through nothing, looked for as such. A circle that touches the
//    first curve at P from one side and grows from nothing first meets each other curve at the
//    radius r(Q) = |Q - P|^2 / (2 (Q - P) . n), least over the points Q of that curve on its
//    side (n the normal at P). Where the two other curves are met at the same radius, the
//    circle touches all three and has none of them inside it. That is one equation in the
//    place of P: its sign changes along the first curve, on either side, are the starts.
//    Where the curves ripple this finds the one circle that sits in the hollow, which the grid
//    of three-point circles does not.
std::vector<Circle> tangentToThree(const Flat &flat) {
    std::vector<Circle> found;
    auto keep = [&](double t[3], Circle circle) {
        if (!settle(flat, t, circle))
            return;
        for (const auto &other : found)
            if (std::hypot(other.x - circle.x, other.y - circle.y) < 1e-7 * circle.radius &&
                std::abs(other.radius - circle.radius) < 1e-7 * circle.radius)
                return;
        found.push_back(circle);
    };
    constexpr int grid = 8;
    for (int i = 0; i < grid; ++i)
        for (int j = 0; j < grid; ++j)
            for (int k = 0; k < grid; ++k) {
                const int index[3] = {i, j, k};
                double t[3];
                At p[3];
                for (int c = 0; c < 3; ++c) {
                    t[c] = flat.first[c] + (flat.last[c] - flat.first[c]) * (index[c] + 0.5) / grid;
                    p[c] = at(flat, c, t[c]);
                }
                Circle circle{};
                if (through(p, circle))
                    keep(t, circle);
            }
    // The second search.
    struct Place {
        double t, x, y;
    };
    std::vector<Place> places[3];
    for (int c = 0; c < 3; ++c) {
        const int n = placesOn(flat.curves[c]);
        for (int i = 0; i <= n; ++i) {
            const double t = flat.first[c] + (flat.last[c] - flat.first[c]) * i / n;
            const gp_Pnt p = flat.curves[c]->Value(t);
            places[c].push_back({t, p.X(), p.Y()});
        }
    }
    for (const double side : {1.0, -1.0}) {
        double last = std::nan("");
        double lastT[3] = {0, 0, 0};
        for (const Place &base : places[0]) {
            const At p = at(flat, 0, base.t);
            const double speed = std::hypot(p.dx, p.dy);
            double value = std::nan(""), t[3] = {base.t, 0, 0}, radius[3] = {0, 0, 0};
            if (speed > 0) {
                const double nx = -p.dy / speed * side, ny = p.dx / speed * side;
                bool both = true;
                for (int c = 1; c < 3; ++c) {
                    double least = std::numeric_limits<double>::infinity();
                    for (const Place &q : places[c]) {
                        const double qx = q.x - p.x, qy = q.y - p.y, across = qx * nx + qy * ny;
                        if (across <= 0)
                            continue;
                        const double r = (qx * qx + qy * qy) / (2 * across);
                        if (r < least) {
                            least = r;
                            t[c] = q.t;
                        }
                    }
                    radius[c] = least;
                    both &= std::isfinite(least);
                }
                if (both)
                    value = radius[1] - radius[2];
            }
            if (!std::isnan(value) && !std::isnan(last) && (value >= 0) != (last >= 0)) {
                // A start from each side of the change: where the nearest point of a curve
                // jumps, only one of them leads anywhere.
                for (const bool here : {true, false}) {
                    double start[3];
                    for (int c = 0; c < 3; ++c)
                        start[c] = here ? t[c] : lastT[c];
                    At q[3];
                    for (int c = 0; c < 3; ++c)
                        q[c] = at(flat, c, start[c]);
                    Circle circle{};
                    if (through(q, circle))
                        keep(start, circle);
                }
            }
            last = value;
            for (int c = 0; c < 3; ++c)
                lastT[c] = t[c];
        }
    }
    return found;
}

// How a circle sits against the curves: whether it is tangent to each of them (its centre at
// its radius from the curve at a place where that distance stops falling or rising, to
// `within`), and whether it cuts through one of them (the curve has points inside the circle
// and points outside it). By OpenCascade's search for those places over each span of each
// curve, not by a fixed number of points.
struct Standing {
    bool tangent, cuts;
};
Standing standing(const Flat &flat, const Circle &circle, double within) {
    const gp_Pnt centre(circle.x, circle.y, 0);
    Standing result{true, false};
    for (size_t c = 0; c < flat.curves.size(); ++c) {
        double least = std::numeric_limits<double>::infinity(), most = 0;
        for (const double t : {flat.first[c], flat.last[c]}) {
            const double d = flat.curves[c]->Value(t).Distance(centre);
            least = std::min(least, d);
            most = std::max(most, d);
        }
        bool tangent = false;
        GeomAPI_ProjectPointOnCurve project(centre, flat.curves[c], flat.first[c], flat.last[c]);
        for (int j = 1; j <= project.NbPoints(); ++j) {
            const double d = project.Distance(j);
            tangent |= std::abs(d - circle.radius) <= within;
            least = std::min(least, d);
            most = std::max(most, d);
        }
        result.tangent &= tangent;
        result.cuts |= least < circle.radius - within && most > circle.radius + within;
    }
    return result;
}

// The curves of the edges in a frame of their own: the plane as z = 0, their middle at the origin.
Flat flatten(const std::vector<Handle(Geom_TrimmedCurve)> &source, const gp_Ax2 &axes) {
    Flat flat;
    gp_Trsf toPlane;
    toPlane.SetTransformation(gp_Ax3(axes));
    Bnd_Box box;
    for (const auto &curve : source)
        BndLib_Add3dCurve::Add(GeomAdaptor_Curve(Handle(Geom_Curve)::DownCast(curve->Transformed(toPlane))), 0, box);
    const gp_Pnt middle((box.CornerMin().XYZ() + box.CornerMax().XYZ()) / 2);
    gp_Trsf toMiddle;
    toMiddle.SetTranslation(gp_Vec(middle.X(), middle.Y(), 0).Reversed());
    flat.toLocal = toMiddle * toPlane;
    flat.extent = box.CornerMin().Distance(box.CornerMax()) / 2;
    for (const auto &curve : source) {
        flat.curves.push_back(Handle(Geom_Curve)::DownCast(curve->Transformed(flat.toLocal)));
        flat.first.push_back(curve->FirstParameter());
        flat.last.push_back(curve->LastParameter());
    }
    return flat;
}
} // namespace

EvolutionData tangentCircle(OcctKernel &kernel, uint32_t first, uint32_t second, uint32_t third,
                            double radius, Vec3 origin, Vec3 normal, Vec3 x, Vec3 near, double tolerance) {
    const char *op = "tangentCircle";
    gp_Ax2 axes(gp_Pnt(origin.x, origin.y, origin.z), gp_Dir(normal.x, normal.y, normal.z),
                gp_Dir(x.x, x.y, x.z));
    const gp_Pln plane(axes);
    std::vector<Handle(Geom_TrimmedCurve)> source;
    auto take = [&](uint32_t id) {
        if (get(kernel, id).ShapeType() != TopAbs_EDGE)
            fail(op, invalid_input, "tangency inputs must be edges");
        double a, b;
        const auto curve = BRep_Tool::Curve(TopoDS::Edge(get(kernel, id)), a, b);
        if (curve.IsNull())
            fail(op, invalid_input, "edge has no 3D curve");
        Handle(Geom_TrimmedCurve) trimmed = new Geom_TrimmedCurve(curve, a, b);
        requirePlanarCurve(op, trimmed, a, b, plane, tolerance);
        source.push_back(trimmed);
    };
    take(first);
    take(second);
    if (third)
        take(third);
    else if (!std::isfinite(radius) || radius <= tolerance)
        fail(op, invalid_input, "positive radius required for two curves");
    const Flat flat = flatten(source, axes);
    const gp_Pnt hint = gp_Pnt(near.x, near.y, near.z).Transformed(flat.toLocal);
    const gp_Pln level(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)));
    auto qualified = [&](int i) {
        return Geom2dGcc_QualifiedCurve(
            Geom2dAdaptor_Curve(GeomAPI::To2d(flat.curves[i], level), flat.first[i], flat.last[i]),
            GccEnt_unqualified);
    };
    // Every circle offered, with how exactly it is to be tangent.
    struct Offer {
        Circle circle;
        double within, distance;
    };
    std::vector<Offer> offers;
    auto offer = [&](const Circle &circle, double within) {
        offers.push_back({circle, within, std::hypot(circle.x - hint.X(), circle.y - hint.Y())});
    };
    if (third) {
        const double exact = acceptedShare * flat.extent;
        Geom2dGcc_Circ2d3Tan solver(qualified(0), qualified(1), qualified(2), tolerance,
                                    (flat.first[0] + flat.last[0]) / 2, (flat.first[1] + flat.last[1]) / 2,
                                    (flat.first[2] + flat.last[2]) / 2);
        if (solver.IsDone())
            for (int i = 1; i <= solver.NbSolutions(); ++i) {
                // OpenCascade's circle is good to its tolerance; from it Newton has the circle
                // itself, also where two of the curves touch and OpenCascade stops short.
                const auto given = solver.ThisSolution(i);
                Circle circle{given.Location().X(), given.Location().Y(), given.Radius()};
                double t[3];
                bool placed = true;
                for (int c = 0; c < 3 && placed; ++c) {
                    GeomAPI_ProjectPointOnCurve project(gp_Pnt(circle.x, circle.y, 0), flat.curves[c],
                                                        flat.first[c], flat.last[c]);
                    double best = std::numeric_limits<double>::infinity();
                    for (int j = 1; j <= project.NbPoints(); ++j)
                        if (std::abs(project.Distance(j) - circle.radius) < best) {
                            best = std::abs(project.Distance(j) - circle.radius);
                            t[c] = project.Parameter(j);
                        }
                    placed = std::isfinite(best);
                }
                Circle settled = circle;
                if (placed && settle(flat, t, settled) &&
                    std::hypot(settled.x - circle.x, settled.y - circle.y) <= 10 * tolerance)
                    offer(settled, exact + exactShare * settled.radius);
                else
                    offer(circle, tolerance);
            }
        for (const auto &circle : tangentToThree(flat))
            offer(circle, exact + exactShare * circle.radius);
    } else {
        Geom2dGcc_Circ2d2TanRad solver(qualified(0), qualified(1), radius, tolerance);
        if (solver.IsDone())
            for (int i = 1; i <= solver.NbSolutions(); ++i) {
                const auto given = solver.ThisSolution(i);
                offer({given.Location().X(), given.Location().Y(), given.Radius()}, tolerance);
            }
    }
    // Of the circles offered, the one whose centre is nearest the hint among those that cut
    // through none of the curves; only if every one cuts, the nearest of those. Each is first
    // checked against the trimmed edges themselves, for inner as for outer contact. They are
    // looked at from the nearest on, so that the search stops at the first that cuts nothing.
    std::sort(offers.begin(), offers.end(), [](const Offer &a, const Offer &b) { return a.distance < b.distance; });
    const Circle *chosen = nullptr;
    for (const auto &candidate : offers) {
        const Standing how = standing(flat, candidate.circle, candidate.within);
        if (!how.tangent)
            continue;
        if (!how.cuts) {
            chosen = &candidate.circle;
            break;
        }
        if (!chosen)
            chosen = &candidate.circle;
    }
    if (!chosen)
        fail(op, no_solution, "no circle tangent to the trimmed curves was found");
    gp_Ax2 location = axes;
    location.SetLocation(gp_Pnt(chosen->x, chosen->y, 0).Transformed(flat.toLocal.Inverted()));
    return stored(kernel,
                  BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(location, chosen->radius)).Edge()).Wire());
}
EvolutionData tangentArcToCurve(OcctKernel &kernel, Vec3 start, Vec3 tangent, uint32_t targetId, Vec3 normal,
                                Vec3 near, double tolerance) {
    const char *op = "tangentArcToCurve";
    const gp_Pnt p(start.x, start.y, start.z);
    const gp_Vec t(tangent.x, tangent.y, tangent.z), n(normal.x, normal.y, normal.z);
    if (t.Magnitude() < 1e-15 || n.Magnitude() < 1e-15 || std::abs(t.Normalized().Dot(n.Normalized())) > BRepExt::tolerance::angle)
        fail(op, invalid_input, "arc tangent must lie in the chosen plane");
    const gp_Ax2 axes(p, gp_Dir(n), gp_Dir(t));
    const gp_Pln plane(axes);
    if (get(kernel, targetId).ShapeType() != TopAbs_EDGE)
        fail(op, invalid_input, "target must be an edge");
    double first, last;
    const auto curve = BRep_Tool::Curve(TopoDS::Edge(get(kernel, targetId)), first, last);
    if (curve.IsNull())
        fail(op, invalid_input, "edge has no 3D curve");
    Handle(Geom_TrimmedCurve) target = new Geom_TrimmedCurve(curve, first, last);
    requirePlanarCurve(op, target, first, last, plane, tolerance);
    const gp_Pnt hint(near.x, near.y, near.z);
    // The arc about `center` from the start, leaving it along the direction, to `q`.
    auto arcTo = [&](const gp_Pnt &center, const gp_Pnt &q, double &angle) {
        const gp_Vec radial(center, p);
        auto axis = gp_Dir(n);
        if (gp_Vec(axis).Crossed(radial).Dot(t) < 0)
            axis.Reverse();
        const gp_Ax2 placement(center, axis, gp_Dir(radial));
        const gp_Vec end(center, q);
        angle = std::atan2(end.Dot(gp_Vec(placement.YDirection())), end.Dot(gp_Vec(placement.XDirection())));
        if (angle <= 0)
            angle += 2 * M_PI;
        return placement;
    };
    // Whether the arc passes through the target before it ends on it. Between two places where
    // the target's distance from the centre stops falling or rising it only falls or only
    // rises, so it is at the radius at most once there: found by halving, and then seen whether
    // that point is on the arc.
    auto crossesBefore = [&](const gp_Pnt &center, double r, const gp_Pnt &q) {
        double angle;
        const gp_Ax2 placement = arcTo(center, q, angle);
        std::vector<double> stops{first, last};
        GeomAPI_ProjectPointOnCurve projection(center, target, first, last);
        for (int j = 1; j <= projection.NbPoints(); ++j)
            stops.push_back(projection.Parameter(j));
        std::sort(stops.begin(), stops.end());
        auto off = [&](double u) { return target->Value(u).Distance(center) - r; };
        for (size_t i = 0; i + 1 < stops.size(); ++i) {
            double lo = stops[i], hi = stops[i + 1];
            const double low = off(lo), high = off(hi);
            if (std::abs(low) <= tolerance || std::abs(high) <= tolerance || (low > 0) == (high > 0))
                continue;
            for (int k = 0; k < 60; ++k) {
                const double mid = (lo + hi) / 2;
                ((off(mid) > 0) == (low > 0) ? lo : hi) = mid;
            }
            const gp_Vec to(center, target->Value((lo + hi) / 2));
            double at = std::atan2(to.Dot(gp_Vec(placement.YDirection())), to.Dot(gp_Vec(placement.XDirection())));
            if (at < 0)
                at += 2 * M_PI;
            if (at * r > tolerance && (angle - at) * r > tolerance)
                return true;
        }
        return false;
    };
    // Of the arcs found, the one that ends nearest the hint among those that do not pass
    // through the target on their way; only if every one does, the nearest of those.
    struct Offer {
        gp_Pnt center, end;
        double radius, distance;
    };
    std::vector<Offer> offers;
    auto offer = [&](const gp_Pnt &center, double r, const gp_Pnt &q) {
        if (r > tolerance && q.Distance(p) > tolerance)
            offers.push_back({center, q, r, q.Distance(hint)});
    };
    Geom2dGcc_QualifiedCurve qualified(Geom2dAdaptor_Curve(GeomAPI::To2d(target, plane)), GccEnt_unqualified);
    Handle(Geom2d_Point) point = new Geom2d_CartesianPoint(0, 0);
    Handle(Geom2d_Curve) centers = new Geom2d_Line(gp_Pnt2d(0, 0), gp_Dir2d(0, 1));
    // OpenCascade's solver has every arc onto a line or a circle; onto a spline it throws.
    try {
        Geom2dGcc_Circ2d2TanOn solver(qualified, point, Geom2dAdaptor_Curve(centers), tolerance,
                                      (first + last) / 2, 0);
        if (solver.IsDone())
            for (int i = 1; i <= solver.NbSolutions(); ++i) {
                const auto c = solver.ThisSolution(i);
                const auto center = p.Translated(gp_Vec(axes.XDirection()) * c.Location().X() +
                                                 gp_Vec(axes.YDirection()) * c.Location().Y());
                GeomAPI_ProjectPointOnCurve projection(center, target, first, last);
                for (int j = 1; j <= projection.NbPoints(); ++j)
                    if (std::abs(projection.Distance(j) - c.Radius()) <= tolerance)
                        offer(center, c.Radius(), projection.Point(j));
            }
    } catch (const Standard_Failure &) {
    }
    // OpenCascade's solver iterates from one start on the target and finds no arc onto most
    // splines. Here every point of the target is tried: the circle that leaves the start in
    // its direction and passes through the point Q(u) has its centre at the signed distance
    // |Q - p|^2 / (2 (Q - p) . m) along m, the normal to the direction; it is tangent to the
    // target where its radius is square to the target there. That is one equation in u: its
    // sign changes over 1024 pieces of the target are closed in on by halving.
    {
        const gp_Vec m = gp_Vec(axes.YDirection());
        auto circleAt = [&](double u, gp_Pnt &center, gp_Pnt &q, double &r) {
            gp_Vec d1;
            target->D1(u, q, d1);
            const gp_Vec to(p, q);
            const double across = to.Dot(m);
            if (std::abs(across) <= 1e-14 * to.Magnitude())
                return std::nan("");
            r = to.SquareMagnitude() / (2 * across);
            center = p.Translated(m * r);
            const double speed = d1.Magnitude();
            return speed < 1e-300 ? std::nan("") : gp_Vec(q, center).Dot(d1) / speed;
        };
        constexpr int pieces = 1024;
        gp_Pnt center, q;
        double r = 0, before = circleAt(first, center, q, r);
        for (int i = 1; i <= pieces; ++i) {
            double lo = first + (last - first) * (i - 1) / pieces, hi = first + (last - first) * i / pieces;
            const double after = circleAt(hi, center, q, r);
            const bool crosses = !std::isnan(before) && !std::isnan(after) && (before >= 0) != (after >= 0);
            const double low = before;
            before = after;
            if (!crosses)
                continue;
            for (int k = 0; k < 60; ++k) {
                const double mid = (lo + hi) / 2, value = circleAt(mid, center, q, r);
                if (std::isnan(value))
                    break;
                ((value >= 0) == (low >= 0) ? lo : hi) = mid;
            }
            const double value = circleAt((lo + hi) / 2, center, q, r);
            // Where the circle through the point turns over (its centre runs off to infinity)
            // the sign changes too, without a tangency: only a true root is taken.
            // A true root: square to the target to 1e-9 of the radius, plus what a double holds
            // of a point this far from the origin.
            if (std::isnan(value) ||
                std::abs(value) > exactShare * std::abs(r) + resolvedShare * (q.XYZ().Modulus() + std::abs(r)))
                continue;
            offer(center, std::abs(r), q);
        }
    }
    std::sort(offers.begin(), offers.end(), [](const Offer &a, const Offer &b) { return a.distance < b.distance; });
    const Offer *chosen = nullptr;
    for (const auto &candidate : offers) {
        if (!crossesBefore(candidate.center, candidate.radius, candidate.end)) {
            chosen = &candidate;
            break;
        }
        if (!chosen)
            chosen = &candidate;
    }
    if (!chosen)
        fail(op, no_solution, "no arc ends tangent to the trimmed target");
    double angle;
    const gp_Ax2 placement = arcTo(chosen->center, chosen->end, angle);
    const double radius = chosen->radius;
    return stored(kernel,
                  BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(placement, radius), 0, angle).Edge()).Wire());
}
} // namespace BRepExt
