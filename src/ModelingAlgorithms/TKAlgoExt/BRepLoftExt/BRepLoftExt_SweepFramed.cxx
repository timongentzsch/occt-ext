// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepLib_FindSurface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <gp_GTrsf.hxx>
#include <gp_Mat.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <vector>

// A sweep with frames of our own. OpenCascade's pipe twists a profile only where the path
// starts square to it and the profile sits on the path, and has no guide that scales.
//
// The frames along the path are rotation-minimising ones, by the double reflection method
// (W. Wang, B. Juettler, D. Zheng, Y. Liu: Computation of rotation minimizing frames, ACM
// Transactions on Graphics 27 (1), 2008): from one station to the next the frame is mirrored
// in the plane between the two points and then in the plane that takes the mirrored tangent
// to the next tangent. Its error is of the fourth order in the step.
//
// The profile is where it was drawn: its station is where the path passes through its plane,
// or, where that is somewhere else along the path or the profile is not flat, the place of the
// path nearest the profile's middle. From there it is carried both ways along the path:
// rigidly with the frame, turned about the path's tangent by the twist, which grows evenly
// with the path's length, and scaled about the path's point by a factor that grows evenly with
// the length and is `scale` times as large at the path's end as at its start. Guides set instead where the section is in each plane square to the path, and twist and
// scale have no part: the point where the first guide crosses the profile's plane stays on it,
// the section otherwise rigid in the path's frame; a second guide turns it about that point,
// and three map the plane affinely (see the guides below).
//
// The result is not skinned through copies of the profile. What carries the profile to a place
// of the path is a map x -> A x + b, and a map that is a cubic spline in the path's parameter
// carries a spline curve to a spline surface exactly: the surface's control points are the
// map's control maps applied to the profile's. So the only thing approximated is the map, as a
// piecewise cubic through its values and slopes at stations, and it is checked where it is
// worst, between the stations: there the spline map may move no corner of the box of the
// profile's control points further than the tolerance from where the true map puts it.
// Stations are doubled until that holds, and until the frames and the path's length have
// stopped changing. The path is taken in pieces over which it has a continuous curvature
// (each gets surfaces of its own, since the map has a kink between them), and each piece in
// stretches between the places where its curvature has a kink (the knots of a cubic spline):
// the slopes at the stations are taken from one stretch only.

namespace loft {
namespace {
const char *op = "sweepFramed";

/// A part of an edge of the path with continuous curvature, run in the path's direction by x
/// from 0 to 1, and the x where the curvature has a kink.
struct Piece {
    BRepAdaptor_Curve curve;
    bool back;
    double first, last;
    std::vector<double> breaks; // rising, from 0 to 1
    double parameter(double x) const { return back ? last - x * (last - first) : first + x * (last - first); }
    double share(double u) const {
        const double x = (u - first) / (last - first);
        return std::clamp(back ? 1 - x : x, 0.0, 1.0);
    }
    void at(double x, gp_Pnt &p, gp_Vec &d) const {
        curve.D1(parameter(x), p, d);
        d *= back ? first - last : last - first;
    }
};
/// A stretch of a piece over which the map is smooth.
struct Stretch {
    int piece;
    double from, to;
};
/// A map of the plane square to the path, x -> A x + d, with A = [[a b] [c e]].
struct PlaneMap {
    double a = 1, b = 0, c = 0, e = 1;
    gp_XY d;
};

struct Node {
    gp_Pnt at;
    gp_Vec t, r;       // tangent and a normal, unit
    double length = 0; // of the path up to here
    PlaneMap guided;   // where the guides put the profile, in the plane of (r, t x r)
};
/// The map at a node: where the three axes of the profile's own frame go, and its origin.
using Map = std::array<gp_Vec, 4>;

/// The normal `r` of the frame (a, ta) carried to (b, tb) by double reflection.
gp_Vec carriedNormal(const gp_Pnt &a, const gp_Vec &ta, const gp_Vec &r, const gp_Pnt &b, const gp_Vec &tb) {
    const gp_Vec v1(a, b);
    const double c1 = v1.SquareMagnitude();
    gp_Vec next = r;
    if (c1 > 0) {
        const gp_Vec rl = r - v1 * (2 / c1 * v1.Dot(r));
        const gp_Vec tl = ta - v1 * (2 / c1 * v1.Dot(ta));
        const gp_Vec v2 = tb - tl;
        const double c2 = v2.SquareMagnitude();
        next = c2 > 0 ? rl - v2 * (2 / c2 * v2.Dot(rl)) : rl;
    }
    // Kept square to the tangent and of length one against rounding.
    next -= tb * next.Dot(tb);
    return next.Normalized();
}

/// The length of a piece between x = a and x = b: Gauss-Legendre with four points.
double lengthOf(const Piece &piece, double a, double b) {
    static const double place[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                    0.8611363115940526};
    static const double weight[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                     0.3478548451374538};
    double sum = 0;
    for (int k = 0; k < 4; ++k) {
        gp_Pnt p;
        gp_Vec d;
        piece.at((a + b) / 2 + (b - a) / 2 * place[k], p, d);
        sum += weight[k] * d.Magnitude();
    }
    return sum * (b - a) / 2;
}

/// The pieces of the path, in its direction.
std::vector<Piece> piecesOf(const TopoDS_Wire &path) {
    std::vector<Piece> pieces;
    for (BRepTools_WireExplorer it(path); it.More(); it.Next()) {
        const BRepAdaptor_Curve curve(it.Current());
        const bool back = it.Current().Orientation() == TopAbs_REVERSED;
        NCollection_Array1<double> smooth(1, curve.NbIntervals(GeomAbs_C2) + 1);
        curve.Intervals(smooth, GeomAbs_C2);
        NCollection_Array1<double> kinks(1, curve.NbIntervals(GeomAbs_C3) + 1);
        curve.Intervals(kinks, GeomAbs_C3);
        std::vector<Piece> ofEdge;
        for (int i = smooth.Lower(); i < smooth.Upper(); ++i) {
            Piece piece{curve, back, smooth(i), smooth(i + 1), {0.0}};
            for (int k = kinks.Lower(); k <= kinks.Upper(); ++k)
                if (kinks(k) > piece.first && kinks(k) < piece.last)
                    piece.breaks.push_back(piece.share(kinks(k)));
            piece.breaks.push_back(1.0);
            std::sort(piece.breaks.begin(), piece.breaks.end());
            ofEdge.push_back(piece);
        }
        if (back)
            std::reverse(ofEdge.begin(), ofEdge.end());
        pieces.insert(pieces.end(), ofEdge.begin(), ofEdge.end());
    }
    return pieces;
}

/// Where the profile is on the path: the piece and x.
struct Station {
    int piece = 0;
    double x = 0;
};
Station stationOf(const TopoDS_Wire &profile, const std::vector<Piece> &pieces, double tolerance) {
    GProp_GProps line;
    BRepGProp::LinearProperties(profile, line);
    const gp_Pnt middle = line.CentreOfMass();
    // How far the profile reaches from its middle.
    double reach = 0;
    for (BRepTools_WireExplorer it(profile); it.More(); it.Next()) {
        Bnd_Box box;
        BRepBndLib::AddOptimal(it.Current(), box, false);
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        for (const double x : {x0, x1})
            for (const double y : {y0, y1})
                for (const double z : {z0, z1})
                    reach = std::max(reach, middle.Distance(gp_Pnt(x, y, z)));
    }
    struct Offer {
        Station at;
        double gap = -1;
    };
    Offer crossing, nearest;
    auto offer = [&](Offer &best, size_t piece, double u, const gp_Pnt &p) {
        const double gap = p.Distance(middle);
        if (best.gap < 0 || gap < best.gap)
            best = {{int(piece), pieces[piece].share(u)}, gap};
    };
    // The profile's plane. Of a copy: the search puts a closed wire into a face of its own,
    // which would mark the caller's wire as part of something.
    BRepLib_FindSurface flat(BRepBuilderAPI_Copy(profile).Shape(), -1, true);
    Handle(Geom_Plane) plane;
    if (flat.Found())
        plane = Handle(Geom_Plane)::DownCast(flat.Surface());
    if (!plane.IsNull())
        plane = new Geom_Plane(plane->Pln().Transformed(flat.Location().Transformation()));
    for (size_t k = 0; k < pieces.size(); ++k) {
        const Piece &piece = pieces[k];
        double first, last;
        const Handle(Geom_Curve) whole = BRep_Tool::Curve(piece.curve.Edge(), first, last);
        if (whole.IsNull())
            continue;
        if (!plane.IsNull()) {
            GeomAPI_IntCS cut(new Geom_TrimmedCurve(whole, piece.first, piece.last), plane);
            if (cut.IsDone())
                for (int i = 1; i <= cut.NbPoints(); ++i) {
                    double u, v, w;
                    cut.Parameters(i, u, v, w);
                    offer(crossing, k, w, cut.Point(i));
                }
        }
        GeomAPI_ProjectPointOnCurve foot(middle, whole, piece.first, piece.last);
        for (int i = 1; i <= foot.NbPoints(); ++i)
            offer(nearest, k, foot.Parameter(i), foot.Point(i));
        for (const double u : {piece.first, piece.last}) {
            // An end of the piece that lies in the plane.
            if (!plane.IsNull() && plane->Pln().Distance(piece.curve.Value(u)) <= tolerance)
                offer(crossing, k, u, piece.curve.Value(u));
            offer(nearest, k, u, piece.curve.Value(u));
        }
    }
    // The profile is where the path passes through its plane, unless that is somewhere else
    // along the path: further from the profile, by more than the profile's size, than the path
    // comes to it anyway. Then, and for a profile that is not flat, it is at the place of the
    // path nearest its middle.
    return crossing.gap >= 0 && crossing.gap <= nearest.gap + reach ? crossing.at : nearest.at;
}

/// Where a guide crosses the plane through `origin` square to `heading`: of the crossings, the
/// one nearest `origin`. An end of the guide within the tolerance of the plane counts.
double crossing(const BRepAdaptor_CompCurve &guide, const gp_Pnt &origin, const gp_Vec &heading,
                double tolerance) {
    const double first = guide.FirstParameter(), last = guide.LastParameter();
    auto height = [&](double w) { return gp_Vec(origin, guide.Value(w)).Dot(heading); };
    double best = -1, at = 0;
    auto offer = [&](double w) {
        const double gap = guide.Value(w).Distance(origin);
        if (best < 0 || gap < best) {
            best = gap;
            at = w;
        }
    };
    for (const double w : {first, last})
        if (std::abs(height(w)) <= tolerance)
            offer(w);
    const int samples = 256;
    double w0 = first, h0 = height(first);
    for (int i = 1; i <= samples; ++i) {
        const double w1 = first + (last - first) * i / samples, h1 = height(w1);
        if (h0 == 0)
            offer(w0);
        else if (h0 * h1 < 0) {
            double lo = w0, hi = w1;
            for (int k = 0; k < 60; ++k) {
                const double mid = (lo + hi) / 2;
                if ((height(mid) < 0) == (h0 < 0))
                    lo = mid;
                else
                    hi = mid;
            }
            offer((lo + hi) / 2);
        }
        w0 = w1;
        h0 = h1;
    }
    if (best < 0)
        fail(op, invalid_input, "a guide does not cross the plane of the profile");
    return at;
}

/// An edge of the profile as a spline, and its control points in the profile's frame.
struct ProfileCurve {
    Handle(Geom_BSplineCurve) curve;
    std::vector<gp_XYZ> local;
};
} // namespace

EvolutionData sweepFramed(OcctKernel &kernel, uint32_t profileId, uint32_t pathId, double twist, double scale,
                          const std::vector<uint32_t> &guideIds, int guideMethod, double tolerance, bool solid) {
    const auto &profile = get(kernel, profileId), &path = get(kernel, pathId);
    if (!std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "tolerance must be positive");
    if (profile.ShapeType() != TopAbs_WIRE || path.ShapeType() != TopAbs_WIRE)
        fail(op, invalid_input, "the profile and the path must be curves");
    if (!std::isfinite(twist) || std::abs(twist) > 8 * M_PI)
        fail(op, invalid_input, "the twist must be within four turns");
    if (!std::isfinite(scale) || scale <= 0)
        fail(op, invalid_input, "the scale must be positive");
    if (guideMethod < 0 || guideMethod > 2)
        fail(op, invalid_input, "the guide method is point (0), chord (1), or curve (2)");
    for (const uint32_t id : guideIds)
        if (get(kernel, id).ShapeType() != TopAbs_WIRE)
            fail(op, invalid_input, "a guide must be a curve");
    if (guideIds.size() > 3)
        fail(op, unsupported_geometry, "at most three guides are followed");
    if (!guideIds.empty() && (twist != 0 || scale != 1))
        fail(op, invalid_input, "the guides set the profile's place, turn and size: no twist or scale with them");
    // Three guides set every point of the plane: the method is for two.
    const int method = guideIds.size() == 2 ? guideMethod : 0;
    const TopoDS_Wire spineWire = TopoDS::Wire(path), profileWire = TopoDS::Wire(profile);
    if (BRep_Tool::IsClosed(spineWire))
        fail(op, unsupported_geometry, "a closed path is not swept with frames of our own");

    const std::vector<Piece> pieces = piecesOf(spineWire);
    if (pieces.empty())
        fail(op, invalid_input, "the path has no edges");
    std::vector<Stretch> stretches;
    std::vector<size_t> firstStretch; // of each piece, and one past the last
    for (size_t k = 0; k < pieces.size(); ++k) {
        firstStretch.push_back(stretches.size());
        for (size_t i = 0; i + 1 < pieces[k].breaks.size(); ++i)
            stretches.push_back({int(k), pieces[k].breaks[i], pieces[k].breaks[i + 1]});
    }
    firstStretch.push_back(stretches.size());

    // The profile's frame: at its station, with a normal carried there below.
    const Station station = stationOf(profileWire, pieces, tolerance);
    const Piece &stationPiece = pieces[size_t(station.piece)];
    gp_Pnt origin;
    gp_Vec heading;
    stationPiece.at(station.x, origin, heading);
    if (heading.Magnitude() <= 0)
        fail(op, invalid_input, "the path has no direction where the profile is");
    heading.Normalize();
    size_t stationStretch = firstStretch[size_t(station.piece)];
    while (stationStretch + 1 < firstStretch[size_t(station.piece) + 1] && stretches[stationStretch].to < station.x)
        ++stationStretch;

    // The profile's edges as splines, and how far it reaches from the path's point.
    std::vector<ProfileCurve> curves;
    double reach = 0;
    for (BRepTools_WireExplorer it(profileWire); it.More(); it.Next()) {
        double first, last;
        const Handle(Geom_Curve) whole = BRep_Tool::Curve(it.Current(), first, last);
        if (whole.IsNull())
            fail(op, invalid_input, "an edge of the profile has no curve");
        Handle(Geom_BSplineCurve) spline =
            GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(whole, first, last));
        if (spline->IsPeriodic())
            spline->SetNotPeriodic();
        for (int i = 1; i <= spline->NbPoles(); ++i)
            reach = std::max(reach, origin.Distance(spline->Pole(i)));
        curves.push_back({spline, {}});
    }
    if (curves.empty() || reach <= 0)
        fail(op, invalid_input, "the profile has no size");
    const bool closedProfile = BRep_Tool::IsClosed(profileWire);

    // The frames turn smoothly only along a path without corners. A kink that moves no point
    // of the profile by the tolerance is none.
    for (size_t k = 0; k + 1 < pieces.size(); ++k) {
        gp_Pnt p;
        gp_Vec out, in;
        pieces[k].at(1, p, out);
        pieces[k + 1].at(0, p, in);
        if (out.Magnitude() <= 0 || in.Magnitude() <= 0)
            fail(op, invalid_input, "the path has no direction somewhere");
        if (out.Angle(in) * reach > tolerance)
            fail(op, unsupported_geometry, "the path has a corner: twist and guides need a smooth path");
    }

    // The guides, and where each crosses the plane of the profile's station.
    std::vector<std::unique_ptr<BRepAdaptor_CompCurve>> guides;
    std::vector<double> crossings;
    for (const uint32_t id : guideIds) {
        const TopoDS_Wire guideWire = TopoDS::Wire(get(kernel, id));
        // At a corner of a guide the section's move has a kink, and the sweep a crease across
        // it: that is not built.
        const std::vector<Piece> guidePieces = piecesOf(guideWire);
        for (size_t k = 0; k + 1 < guidePieces.size(); ++k) {
            gp_Pnt p;
            gp_Vec out, in;
            guidePieces[k].at(1, p, out);
            guidePieces[k + 1].at(0, p, in);
            if (out.Magnitude() > 0 && in.Magnitude() > 0 && out.Angle(in) * reach > tolerance)
                fail(op, unsupported_geometry, "a guide has a corner: a guide has to be a smooth curve");
        }
        guides.push_back(std::make_unique<BRepAdaptor_CompCurve>(guideWire, false));
        crossings.push_back(crossing(*guides.back(), origin, heading, tolerance));
        gp_Pnt g;
        gp_Vec dg;
        guides.back()->D1(crossings.back(), g, dg);
        if (std::abs(dg.Dot(heading)) <= loft::tolerance::angle * dg.Magnitude())
            fail(op, invalid_input, "a guide runs across the path where the profile is, not beside it");
    }

    const double twistWeight = std::abs(twist) * reach; // how far a share of the twist moves a point
    std::vector<TopoDS_Face> faces;
    std::array<std::pair<gp_Mat, gp_XYZ>, 2> ends; // the map x -> A x + b of each end
    gp_Vec endNormalBefore;
    double lengthBefore = -1;
    for (int stations = 4;; stations *= 2) {
        checkCancelled(op);
        if (stations > 2048)
            fail(op, no_solution,
                 "the sweep could not be made to tolerance with 2048 stations to each stretch of the path");
        const int fine = 4 * stations; // steps to a stretch; the stations are every fourth node
        const auto count = size_t(fine) + 1;
        std::vector<Node> nodes(stretches.size() * count);
        auto xOf = [&](size_t s, int j) {
            return stretches[s].from + (stretches[s].to - stretches[s].from) * j / fine;
        };
        double length = 0;
        for (size_t s = 0; s < stretches.size(); ++s)
            for (int j = 0; j <= fine; ++j) {
                Node &node = nodes[s * count + size_t(j)];
                const Piece &piece = pieces[size_t(stretches[s].piece)];
                gp_Vec d;
                piece.at(xOf(s, j), node.at, d);
                if (d.Magnitude() <= 0)
                    fail(op, invalid_input, "the path has no direction somewhere");
                node.t = d.Normalized();
                if (j > 0)
                    length += lengthOf(piece, xOf(s, j - 1), xOf(s, j));
                node.length = length;
            }
        if (length <= tolerance)
            fail(op, invalid_input, "the path has no length");
        {
            const gp_Vec t0 = nodes[0].t;
            nodes[0].r = t0.Crossed(std::abs(t0.X()) < 0.9 ? gp_Vec(1, 0, 0) : gp_Vec(0, 1, 0)).Normalized();
        }
        for (size_t i = 1; i < nodes.size(); ++i)
            nodes[i].r = carriedNormal(nodes[i - 1].at, nodes[i - 1].t, nodes[i - 1].r, nodes[i].at, nodes[i].t);
        // Have the frames and the length stopped changing? (The first normal is the same each time.)
        const gp_Vec endNormal = nodes.back().r;
        const bool settled =
            lengthBefore >= 0 &&
            (endNormalBefore - endNormal).Magnitude() * reach * std::max(1.0, scale) <= tolerance / 4 &&
            std::abs(length - lengthBefore) * (twistWeight + std::abs(scale - 1) * reach) <= tolerance / 4 * length;
        endNormalBefore = endNormal;
        lengthBefore = length;

        // The profile's frame, carried from the node before its station.
        const Stretch &home = stretches[stationStretch];
        const int homeStep =
            std::clamp(int(std::floor((station.x - home.from) / (home.to - home.from) * fine)), 0, fine);
        const size_t before = stationStretch * count + size_t(homeStep);
        const double stationLength =
            nodes[before].length + lengthOf(stationPiece, xOf(stationStretch, homeStep), station.x);
        const gp_Vec normal = carriedNormal(nodes[before].at, nodes[before].t, nodes[before].r, origin, heading);
        const gp_Vec side = heading.Crossed(normal);
        auto local = [&](const gp_Pnt &p) {
            const gp_Vec v(origin, p);
            return gp_XYZ(v.Dot(normal), v.Dot(side), v.Dot(heading));
        };
        const auto sizeAt = [&](double along) { return 1 + (scale - 1) * along / length; };

        // The profile's control points in its frame, and their box.
        gp_XYZ low(1e300, 1e300, 1e300), high(-1e300, -1e300, -1e300);
        for (auto &c : curves) {
            c.local.clear();
            for (int i = 1; i <= c.curve->NbPoles(); ++i) {
                const gp_XYZ x = local(c.curve->Pole(i));
                c.local.push_back(x);
                for (int a = 1; a <= 3; ++a) {
                    low.SetCoord(a, std::min(low.Coord(a), x.Coord(a)));
                    high.SetCoord(a, std::max(high.Coord(a), x.Coord(a)));
                }
            }
        }
        // Where the guides put the profile at each node, marched both ways from the profile. The
        // first guide's point keeps its place on that guide, and the profile its size and the
        // path's frame. A second guide turns the section about that point too, so that the
        // second guide's point stays on it and the section is scaled to keep it there (Point),
        // points at it (Chord), or so that the profile passes through it, the point of contact
        // sliding along the profile (Curve). Three set the section by the map of the plane that
        // takes their three points to where those guides are.
        if (!guides.empty()) {
            std::vector<gp_XY> contacts;
            for (size_t k = 0; k < guides.size(); ++k) {
                const gp_XYZ x = local(guides[k]->Value(crossings[k]));
                contacts.emplace_back(x.X(), x.Y());
            }
            const gp_XY &pin = contacts[0];
            if (guides.size() == 2 && (contacts[1] - pin).Modulus() <= tolerance)
                fail(op, invalid_input, "the two guides cross the profile's plane at one point");
            const gp_XY spanA = guides.size() == 3 ? contacts[1] - pin : gp_XY(),
                        spanB = guides.size() == 3 ? contacts[2] - pin : gp_XY();
            const double det = spanA.Crossed(spanB);
            if (guides.size() == 3 && std::abs(det) <= tolerance * std::max(spanA.Modulus(), spanB.Modulus()))
                fail(op, invalid_input, "the three guides cross the profile's plane in a line");
            if (method == 2)
                for (const auto &c : curves)
                    for (const gp_XYZ &x : c.local)
                        if (std::abs(x.Z()) > tolerance)
                            fail(op, unsupported_geometry, "Curve with two guides needs a flat profile square to the path");
            // For Curve: the directions from the pin of the profile's points at `gap` from it.
            auto directionsAt = [&](double gap) {
                std::vector<double> found;
                for (const auto &c : curves) {
                    const double first = c.curve->FirstParameter(), last = c.curve->LastParameter();
                    auto at = [&](double u) {
                        const gp_XYZ x = local(c.curve->Value(u));
                        return gp_XY(x.X(), x.Y()) - pin;
                    };
                    auto off = [&](double u) { return at(u).Modulus() - gap; };
                    auto keep = [&](double u) { found.push_back(std::atan2(at(u).Y(), at(u).X())); };
                    const int samples = 32 * c.curve->NbPoles();
                    std::vector<double> u(size_t(samples) + 1), f(size_t(samples) + 1);
                    for (int i = 0; i <= samples; ++i) {
                        u[size_t(i)] = first + (last - first) * i / samples;
                        f[size_t(i)] = off(u[size_t(i)]);
                    }
                    for (size_t i = 0; i + 1 < u.size(); ++i)
                        if (f[i] == 0)
                            keep(u[i]);
                        else if (f[i] * f[i + 1] < 0) {
                            double lo = u[i], hi = u[i + 1];
                            for (int k = 0; k < 60; ++k) {
                                const double mid = (lo + hi) / 2;
                                if ((off(mid) < 0) == (f[i] < 0))
                                    lo = mid;
                                else
                                    hi = mid;
                            }
                            keep((lo + hi) / 2);
                        }
                    // Where the profile only comes to `gap` and turns back, there is no change of
                    // sign: the nearest it comes, if within the tolerance.
                    for (size_t i = 1; i + 1 < u.size(); ++i)
                        if (std::abs(f[i]) <= std::abs(f[i - 1]) && std::abs(f[i]) <= std::abs(f[i + 1]) &&
                            f[i - 1] * f[i + 1] > 0) {
                            double lo = u[i - 1], hi = u[i + 1];
                            for (int k = 0; k < 60; ++k) {
                                const double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
                                if (std::abs(off(m1)) < std::abs(off(m2)))
                                    hi = m2;
                                else
                                    lo = m1;
                            }
                            if (std::abs(off((lo + hi) / 2)) <= tolerance)
                                keep((lo + hi) / 2);
                        }
                }
                return found;
            };
            // The turn about the pin for Curve nearest `previous`, that puts the profile through `to`.
            auto curveTurn = [&](const gp_XY &to, double previous) {
                double best = 0, gap = -1;
                for (const double direction : directionsAt(to.Modulus())) {
                    const double turn = nearPeriod(std::atan2(to.Y(), to.X()) - direction, previous, 2 * M_PI);
                    if (gap < 0 || std::abs(turn - previous) < gap) {
                        gap = std::abs(turn - previous);
                        best = turn;
                    }
                }
                if (gap < 0)
                    fail(op, no_solution, "the profile does not reach the second guide: Curve turns it but keeps its size");
                return best;
            };
            if (method == 2 && std::abs(curveTurn(contacts[1] - pin, 0)) * (contacts[1] - pin).Modulus() > tolerance)
                fail(op, invalid_input, "for Curve the second guide has to touch the profile");
            auto mapFor = [&](const std::vector<gp_XY> &to, double &turn) {
                PlaneMap m;
                if (guides.size() == 2) {
                    const gp_XY from = contacts[1] - pin, toward = to[1] - to[0];
                    turn = method == 2 ? curveTurn(toward, turn)
                                       : nearPeriod(std::atan2(toward.Y(), toward.X()) - std::atan2(from.Y(), from.X()),
                                                    turn, 2 * M_PI);
                    const double factor = method == 0 ? toward.Modulus() / from.Modulus() : 1;
                    m.a = m.e = factor * std::cos(turn);
                    m.c = factor * std::sin(turn);
                    m.b = -m.c;
                } else if (guides.size() == 3) {
                    // The map that takes the spans from the pin to where the guides are.
                    const gp_XY toA = to[1] - to[0], toB = to[2] - to[0];
                    m.a = (toA.X() * spanB.Y() - toB.X() * spanA.Y()) / det;
                    m.b = (toB.X() * spanA.X() - toA.X() * spanB.X()) / det;
                    m.c = (toA.Y() * spanB.Y() - toB.Y() * spanA.Y()) / det;
                    m.e = (toB.Y() * spanA.X() - toA.Y() * spanB.X()) / det;
                }
                m.d = to[0] - gp_XY(m.a * pin.X() + m.b * pin.Y(), m.c * pin.X() + m.e * pin.Y());
                return m;
            };
            std::vector<double> senses;
            for (size_t k = 0; k < guides.size(); ++k) {
                gp_Pnt g;
                gp_Vec dg;
                guides[k]->D1(crossings[k], g, dg);
                senses.push_back(dg.Dot(heading));
            }
            for (const int step : {1, -1}) {
                std::vector<double> w = crossings;
                double turn = 0;
                for (auto i = std::ptrdiff_t(before) + (step > 0 ? 1 : 0);
                     i >= 0 && i < std::ptrdiff_t(nodes.size()); i += step) {
                    Node &node = nodes[size_t(i)];
                    const gp_Vec binormal = node.t.Crossed(node.r);
                    std::vector<gp_XY> to;
                    for (size_t k = 0; k < guides.size(); ++k) {
                        const BRepAdaptor_CompCurve &guide = *guides[k];
                        // The guide's point in the plane square to the path here.
                        gp_Pnt g;
                        gp_Vec dg;
                        bool found = false;
                        for (int it = 0; it < 50 && !found; ++it) {
                            guide.D1(w[k], g, dg);
                            const double h = gp_Vec(node.at, g).Dot(node.t), slope = dg.Dot(node.t);
                            if (slope * senses[k] <= 0)
                                fail(op, no_solution, "a guide doubles back along the path: it has to run beside it one way");
                            const double move = -h / slope;
                            if (std::abs(move) * dg.Magnitude() <= 1e-3 * tolerance) {
                                found = true;
                                break;
                            }
                            const double next = std::clamp(w[k] + move, guide.FirstParameter(), guide.LastParameter());
                            // Already at the guide's end and asked to go on: it is too short.
                            if (next == w[k])
                                fail(op, no_solution, "a guide does not run beside the whole path");
                            w[k] = next;
                        }
                        if (!found)
                            fail(op, no_solution, "a guide could not be followed along the path");
                        const gp_Vec d(node.at, g);
                        to.emplace_back(d.Dot(node.r), d.Dot(binormal));
                    }
                    node.guided = mapFor(to, turn);
                }
            }
        }

        // The map at each node.
        std::vector<Map> maps(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            const Node &node = nodes[i];
            const gp_Vec binormal = node.t.Crossed(node.r);
            if (!guides.empty()) {
                const PlaneMap &m = node.guided;
                maps[i] = {node.r * m.a + binormal * m.c, node.r * m.b + binormal * m.e, node.t,
                           gp_Vec(node.at.XYZ()) + node.r * m.d.X() + binormal * m.d.Y()};
                continue;
            }
            const double angle = twist * (node.length - stationLength) / length,
                         factor = sizeAt(node.length) / sizeAt(stationLength);
            maps[i] = {(node.r * std::cos(angle) + binormal * std::sin(angle)) * factor,
                       (binormal * std::cos(angle) - node.r * std::sin(angle)) * factor, node.t * factor,
                       gp_Vec(node.at.XYZ())};
        }

        auto applied = [](const Map &m, const gp_XYZ &x) { return m[3] + m[0] * x.X() + m[1] * x.Y() + m[2] * x.Z(); };

        // Slopes by x at every node, by differences of the fourth order within a stretch;
        // where two stretches of a piece meet the map has one slope, their mean.
        std::vector<Map> slopes(nodes.size());
        for (size_t s = 0; s < stretches.size(); ++s) {
            const Map *m = &maps[s * count];
            const double h = (stretches[s].to - stretches[s].from) / fine;
            for (int j = 0; j <= fine; ++j)
                for (size_t a = 0; a < 4; ++a) {
                    gp_Vec d;
                    if (j >= 2 && j <= fine - 2)
                        d = (m[j - 2][a] - m[j - 1][a] * 8 + m[j + 1][a] * 8 - m[j + 2][a]) / (12 * h);
                    else if (j < 2)
                        d = (m[j][a] * -25 + m[j + 1][a] * 48 - m[j + 2][a] * 36 + m[j + 3][a] * 16 - m[j + 4][a] * 3) /
                            (12 * h);
                    else
                        d = (m[j][a] * 25 - m[j - 1][a] * 48 + m[j - 2][a] * 36 - m[j - 3][a] * 16 + m[j - 4][a] * 3) /
                            (12 * h);
                    slopes[s * count + size_t(j)][a] = d;
                }
            if (s > 0 && stretches[s - 1].piece == stretches[s].piece)
                for (size_t a = 0; a < 4; ++a) {
                    const gp_Vec mean = (slopes[s * count][a] + slopes[s * count - 1][a]) / 2;
                    slopes[s * count][a] = slopes[s * count - 1][a] = mean;
                }
        }
        // The cubic through values and slopes at the stations, against the map at the nodes
        // between them; and that the profile moves forwards everywhere.
        bool near = true, forward = true;
        for (size_t s = 0; s < stretches.size(); ++s) {
            const Map *m = &maps[s * count];
            const double span = (stretches[s].to - stretches[s].from) / stations;
            for (int j = 0; j <= fine; ++j) {
                const size_t i = s * count + size_t(j);
                if (j % 4 != 0) {
                    const int j0 = j - j % 4;
                    const double u = double(j % 4) / 4;
                    const Map &y0 = m[j0], &y1 = m[j0 + 4], &d0 = slopes[s * count + size_t(j0)],
                              &d1 = slopes[s * count + size_t(j0) + 4];
                    Map off;
                    for (size_t a = 0; a < 4; ++a)
                        off[a] = y0[a] * (2 * u * u * u - 3 * u * u + 1) + d0[a] * (span * (u * u * u - 2 * u * u + u)) +
                                 y1[a] * (-2 * u * u * u + 3 * u * u) + d1[a] * (span * (u * u * u - u * u)) - m[j][a];
                    for (const double x : {low.X(), high.X()})
                        for (const double y : {low.Y(), high.Y()})
                            for (const double z : {low.Z(), high.Z()})
                                if (applied(off, gp_XYZ(x, y, z)).Magnitude() > tolerance)
                                    near = false;
                }
                // Every point of the profile has to move forwards along the path: where one
                // stands still or moves back, the surface folds over itself there.
                const Map &d = slopes[i];
                const gp_Vec &t = nodes[i].t;
                const double base = d[3].Dot(t);
                const gp_XYZ lean(d[0].Dot(t), d[1].Dot(t), d[2].Dot(t));
                for (const auto &c : curves) {
                    bool hull = true;
                    for (const auto &x : c.local)
                        if (base + lean.Dot(x) <= 0)
                            hull = false;
                    if (hull)
                        continue;
                    // The control points say only that it may: ask the curve itself.
                    const int samples = 32 * c.curve->NbPoles();
                    for (int k = 0; k <= samples && forward; ++k) {
                        const double w = c.curve->FirstParameter() +
                                         (c.curve->LastParameter() - c.curve->FirstParameter()) * k / samples;
                        if (base + lean.Dot(local(c.curve->Value(w))) <= 0)
                            forward = false;
                    }
                }
            }
        }
        if (!forward)
            fail(op, no_solution,
                 "the sweep folds over itself: the profile reaches past the centre of the path's bend, or the "
                 "guide turns or shrinks it faster than it moves on");
        if (!near || !settled)
            continue;

        // The surfaces: for each piece of the path and each edge of the profile, the control
        // points of the map's cubic applied to the edge's control points.
        faces.clear();
        for (size_t k = 0; k < pieces.size(); ++k) {
            const size_t from = firstStretch[k], to = firstStretch[k + 1];
            const int spans = int(to - from) * stations;
            std::vector<Map> control;
            control.reserve(size_t(2 * spans + 2));
            NCollection_Array1<double> vKnots(1, spans + 1);
            NCollection_Array1<int> vMults(1, spans + 1);
            int knot = 0;
            for (size_t s = from; s < to; ++s) {
                const double span = (stretches[s].to - stretches[s].from) / stations;
                for (int j = 0; j <= fine; j += 4) {
                    const Map &y = maps[s * count + size_t(j)], &d = slopes[s * count + size_t(j)];
                    Map in, out;
                    for (size_t a = 0; a < 4; ++a) {
                        in[a] = y[a] - d[a] * (span / 3);
                        out[a] = y[a] + d[a] * (span / 3);
                    }
                    // Into a station from the stretch before it, out of it into the one after.
                    if (j > 0)
                        control.push_back(in);
                    else if (s == from)
                        control.push_back(y);
                    if (j < fine)
                        control.push_back(out);
                    else if (s + 1 == to)
                        control.push_back(y);
                    if (j > 0 || s == from) {
                        ++knot;
                        vKnots(knot) = double(k) + xOf(s, j);
                        vMults(knot) = 2;
                    }
                }
            }
            vMults(1) = vMults(spans + 1) = 4;
            const int rows = int(control.size());
            for (const auto &c : curves) {
                const int columns = c.curve->NbPoles();
                NCollection_Array2<gp_Pnt> poles(1, columns, 1, rows);
                NCollection_Array2<double> weights(1, columns, 1, rows);
                for (int i = 1; i <= columns; ++i)
                    for (int j = 1; j <= rows; ++j) {
                        poles(i, j) = gp_Pnt(applied(control[size_t(j - 1)], c.local[size_t(i - 1)]).XYZ());
                        weights(i, j) = c.curve->Weight(i);
                    }
                const Handle(Geom_BSplineSurface) surface =
                    new Geom_BSplineSurface(poles, weights, c.curve->Knots(), vKnots, c.curve->Multiplicities(),
                                            vMults, c.curve->Degree(), 3);
                BRepBuilderAPI_MakeFace face(surface, tolerance);
                if (!face.IsDone())
                    fail(op, no_solution, "a face of the sweep could not be made");
                faces.push_back(face.Face());
            }
        }
        // Where the ends are, for the caps: the maps there, taken from the profile's frame to the
        // world's.
        for (const bool last : {false, true}) {
            const Map &map = last ? maps.back() : maps.front();
            gp_Mat linear;
            for (int i = 1; i <= 3; ++i)
                for (int j = 1; j <= 3; ++j)
                    linear.SetValue(i, j, map[0].Coord(i) * normal.Coord(j) + map[1].Coord(i) * side.Coord(j) +
                                              map[2].Coord(i) * heading.Coord(j));
            ends[last] = {linear, map[3].XYZ() - origin.XYZ().Multiplied(linear)};
        }
        break;
    }

    // Flat ends close a closed profile into a solid.
    bool capped = false;
    if (solid && closedProfile) {
        std::vector<TopoDS_Face> caps;
        for (const auto &[linear, shift] : ends) {
            // A copy first: a profile that is carried nowhere would be put into the cap as it
            // is, which marks the caller's own wire as part of something.
            const TopoDS_Shape copy = BRepBuilderAPI_Copy(profile).Shape();
            // Turned and scaled alike in every direction, the profile's curves keep their kind;
            // three guides may stretch it more one way than another.
            const gp_Mat square = linear.Transposed() * linear;
            const double size = square.Value(1, 1);
            bool similar = true;
            for (int i = 1; i <= 3; ++i)
                for (int j = 1; j <= 3; ++j)
                    if (std::abs(square.Value(i, j) - (i == j ? size : 0)) > 1e-12 * size)
                        similar = false;
            TopoDS_Shape placed;
            if (similar) {
                gp_Trsf place;
                place.SetValues(linear(1, 1), linear(1, 2), linear(1, 3), shift.X(), linear(2, 1), linear(2, 2),
                                linear(2, 3), shift.Y(), linear(3, 1), linear(3, 2), linear(3, 3), shift.Z());
                placed = BRepBuilderAPI_Transform(copy, place, false).Shape();
            } else
                placed = BRepBuilderAPI_GTransform(copy, gp_GTrsf(linear, shift), false).Shape();
            BRepBuilderAPI_MakeFace cap(TopoDS::Wire(placed), true);
            if (cap.IsDone())
                caps.push_back(cap.Face());
        }
        if (caps.size() == 2) {
            faces.insert(faces.end(), caps.begin(), caps.end());
            capped = true;
        }
    }
    BRepBuilderAPI_Sewing sewing(tolerance);
    for (const auto &face : faces)
        sewing.Add(face);
    sewing.Perform();
    TopoDS_Shape result = sewing.SewedShape();
    if (result.IsNull())
        fail(op, no_solution, "the sweep's faces could not be joined");
    if (capped) {
        if (sewing.NbFreeEdges() > 0)
            fail(op, no_solution, "the sweep's faces did not close into a solid");
        TopExp_Explorer shells(result, TopAbs_SHELL);
        if (!shells.More())
            fail(op, no_solution, "the sweep's faces did not close into a solid");
        BRepBuilderAPI_MakeSolid solid(TopoDS::Shell(shells.Current()));
        shells.Next();
        if (shells.More() || !solid.IsDone())
            fail(op, no_solution, "the sweep's faces did not close into one solid");
        TopoDS_Solid body = solid.Solid();
        BRepLib::OrientClosedSolid(body);
        result = body;
    }
    // A check of the topology only: that the sweep does not fold is established above, and
    // that it does not pass through itself somewhere else along the path is not claimed.
    if (!BRepCheck_Analyzer(result).IsValid())
        fail(op, no_solution, "the sweep's faces do not make a valid body");
    return stored(kernel, result);
}
} // namespace loft
