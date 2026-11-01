// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Whole-operation comparison of TKAlgoExt against the closest stock OCCT facility, on inputs
// with a closed-form answer. CPU time is std::clock() (process CPU), one warm-up and ten
// measured runs per row. Build and run with adm/scripts/extensions/operation-benchmark.sh.
//
// Rows (each stock call is made with the same input as the extension call):
//   fillet  a plate with one through hole, the hole's top edge rounded: BRepFilletAPI_MakeFillet
//           against filletRange. Quality: the volume removed has the closed form
//           follows from Pappus's theorem, see removedByHoleFillet.
//   patch   four cubic Bezier sides of one analytic surface: BRepFill_Filling (G0 on the sides)
//           against patch. Quality: the largest distance of the result from the analytic surface
//           and from the sides, sampled on 33 x 33 and 129 points.
//   hollow  a box with one rounded edge, hollowed to a wall t: BRepOffsetAPI_MakeThickSolid
//           against hollow. Quality: validity (BRepCheck_Analyzer) and the closed-form volume.
//   offset  a rectangle's outline moved outwards by d with arcs at the corners:
//           BRepOffsetAPI_MakeOffset against offsetCurve. Quality: the length P + 2 pi d.

#include <BRepToolsExt_Operations.hxx>
#include <BRepBlendExt_Request.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <BRepPatchExt_Side.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomConvert.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_List.hxx>
#include <GeomLib.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <functional>
#include <string>
#include <typeinfo>
#include <vector>

// A host that keeps shapes in a vector, as BRepToolsExt_Session.hxx asks.
namespace {
std::vector<TopoDS_Shape> shapes;
alignas(8) char token[8];
OcctKernel& session() { return *reinterpret_cast<OcctKernel*>(token); }
} // namespace
namespace BRepExt {
uint32_t store(OcctKernel&, const TopoDS_Shape& shape) {
  shapes.push_back(shape);
  return static_cast<uint32_t>(shapes.size());
}
const TopoDS_Shape& get(const OcctKernel&, uint32_t id) {
  if (id == 0 || id > shapes.size()) fail("shape", invalid_input, "Invalid shape reference");
  return shapes[id - 1];
}
void release(OcctKernel&, uint32_t) {}
double booleanNow() { return 1000. * std::clock() / CLOCKS_PER_SEC; }
bool booleanTimingEnabled() { return false; }
void booleanPhase(const char*, double) {}
namespace patch_host {
void chain(std::vector<Side>&, double) { fail("patch", unsupported_geometry, "open boundaries need a host"); }
void closeGaps(const std::vector<Side>&, std::vector<Rim>&) {}
} // namespace patch_host
} // namespace BRepExt

namespace {

// Checking a result (validity, volume, length) is not part of the operation: those functions stop the
// clock for as long as they run.
std::clock_t pausedTicks = 0;
struct Pause {
  std::clock_t begin = std::clock();
  ~Pause() { pausedTicks += std::clock() - begin; }
};

struct Stats {
  double median = 0, low = 0, high = 0;
  int ok = 0, runs = 0;
  std::string refusal; // the last refusal's text, for a row with no success
};

// Run `work` once to warm up and ten times for the statistics; `work` returns false on a refusal.
Stats measure(const std::function<bool()>& work) {
  Stats s;
  std::vector<double> times;
  for (int i = 0; i < 11; ++i) {
    const std::clock_t begin = std::clock();
    pausedTicks = 0;
    bool ok = false;
    try { ok = work(); }
    catch (const Standard_Failure& e) { ok = false; s.refusal = std::string(typeid(e).name()) + ": " + (e.GetMessageString() ? e.GetMessageString() : ""); }
    catch (const std::exception& e) { ok = false; s.refusal = std::string("std::exception: ") + e.what(); }
    catch (...) { ok = false; s.refusal = "an exception of unknown type"; }
    const double ms = 1000. * (std::clock() - begin - pausedTicks) / CLOCKS_PER_SEC;
    if (i == 0) continue;
    ++s.runs;
    if (ok) { ++s.ok; times.push_back(ms); }
  }
  if (!times.empty()) {
    std::sort(times.begin(), times.end());
    s.low = times.front();
    s.high = times.back();
    s.median = times[times.size() / 2];
  }
  return s;
}

// The volume by adaptive Gauss-Kronrod integration to a relative accuracy of 1e-8. The default
// BRepGProp::VolumeProperties uses a fixed rule, which is off by about 1e-5 of the volume on spline
// blend faces and would be mistaken for a difference between two implementations.
double volume(const TopoDS_Shape& shape) {
  Pause pause;
  GProp_GProps props;
  BRepGProp::VolumePropertiesGK(shape, props, 1e-8, false, true);
  return props.Mass();
}
double length(const TopoDS_Shape& shape) {
  Pause pause;
  GProp_GProps props;
  BRepGProp::LinearProperties(shape, props);
  return props.Mass();
}
bool valid(const TopoDS_Shape& shape) { Pause pause; return !shape.IsNull() && BRepCheck_Analyzer(shape).IsValid(); }
// The same check, left on the clock: a call that must return only valid bodies pays for it.
bool validTimed(const TopoDS_Shape& shape) { return !shape.IsNull() && BRepCheck_Analyzer(shape).IsValid(); }

void row(const char* operation, const char* input, const char* who, const Stats& s, const char* quality) {
  const std::string refused = "REFUSED: " + s.refusal;
  if (s.ok == 0) quality = refused.c_str();
  std::printf("%-8s | %-34s | %-9s | %2d/%2d | %9.3f [%9.3f, %9.3f] | %s\n", operation, input, who, s.ok, s.runs,
              s.median, s.low, s.high, quality);
}

// ---- fillet ---------------------------------------------------------------------------------
// A plate 100 x 80 x 2 with a through hole of radius 1 about its centre; the hole's top edge is
// rounded by r. In the half-plane (rho, z) about the hole axis the rounding removes the square
// [1, 1 + r] x [2 - r, 2] minus the quarter disc of radius r centred at (1 + r, 2 - r). Its area
// is r^2 (1 - pi/4); the square has centroid rho = 1 + r/2 and the quarter disc rho = 1 + r -
// 4 r / (3 pi). Pappus (volume = 2 pi * centroid radius * area) gives the volume removed.
double removedByHoleFillet(double r) {
  const double square = r * r, disc = M_PI * r * r / 4;
  const double squareRho = 1 + r / 2, discRho = 1 + r - 4 * r / (3 * M_PI);
  return 2 * M_PI * (squareRho * square - discRho * disc);
}

void fillet(double r) {
  const double plate = 100 * 80 * 2 - M_PI * 1 * 1 * 2;
  TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Pnt(-50, -40, 0), 100, 80, 2).Shape();
  TopoDS_Shape hole = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -1), gp_Dir(0, 0, 1)), 1, 4).Shape();
  const TopoDS_Shape body = BRepAlgoAPI_Cut(box, hole).Shape();
  TopoDS_Edge top;
  for (TopExp_Explorer it(body, TopAbs_EDGE); it.More(); it.Next()) {
    const TopoDS_Edge edge = TopoDS::Edge(it.Current());
    BRepAdaptor_Curve curve(edge);
    if (curve.GetType() == GeomAbs_Circle && std::abs(curve.Circle().Location().Z() - 2) < 1e-9) top = edge;
  }
  const double expected = plate - removedByHoleFillet(r);
  char input[64];
  std::snprintf(input, sizeof input, "plate, hole r=1, round r=%.2f", r);
  double errorStock = -1, errorExt = -1, volumeStock = 0, volumeExt = 0;
  const Stats stock = measure([&] {
    BRepFilletAPI_MakeFillet maker(body);
    maker.Add(r, top);
    maker.Build();
    if (!maker.IsDone() || !valid(maker.Shape())) return false;
    volumeStock = volume(maker.Shape());
    errorStock = std::abs(volumeStock - expected);
    return true;
  });
  const Stats checked = measure([&] {
    BRepFilletAPI_MakeFillet maker(body);
    maker.Add(r, top);
    maker.Build();
    return maker.IsDone() && validTimed(maker.Shape());
  });
  const uint32_t bodyId = BRepExt::store(session(), body), edgeId = BRepExt::store(session(), top);
  const Stats ext = measure([&] {
    const auto result = BRepExt::filletRange(session(), bodyId, edgeId, r, 0, 1, true, {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    if (!valid(made)) return false;
    volumeExt = volume(made);
    errorExt = std::abs(volumeExt - expected);
    return true;
  });
  char q[160];
  std::snprintf(q, sizeof q, "V %.12f (expected %.12f), |dV| %.3e", volumeStock, expected, errorStock);
  row("fillet", input, "stock", stock, q);
  row("fillet", input, "stock+check", checked, q);
  std::snprintf(q, sizeof q, "V %.12f (expected %.12f), |dV| %.3e", volumeExt, expected, errorExt);
  row("fillet", input, "extension", ext, q);
}

// ---- patch ----------------------------------------------------------------------------------
// The surface S(u,v) = (u, v, a (u^3 - 3 u v^2)) on [0,1]^2 (a harmonic cubic) with its four
// sides as cubic Bezier edges; the cubic is exactly its own Bezier curve on each side.
gp_Pnt surfacePoint(double u, double v, double a) { return gp_Pnt(u, v, a * (u * u * u - 3 * u * v * v)); }

Handle(Geom_BezierCurve) side(double a, bool alongU, double fixed) {
  // Exact cubic Bezier through the monomial data of the side: sample 4 Bernstein coefficients by
  // interpolating at t = 0, 1/3, 2/3, 1 (a cubic is reproduced exactly).
  const auto point = [&](double t) { return alongU ? surfacePoint(t, fixed, a) : surfacePoint(fixed, t, a); };
  const gp_Pnt p0 = point(0), p1 = point(1. / 3), p2 = point(2. / 3), p3 = point(1);
  // Solve for the Bezier poles B1, B2 from the interpolation conditions at 1/3 and 2/3.
  const gp_XYZ b1 = (-5 * p0.XYZ() + 18 * p1.XYZ() - 9 * p2.XYZ() + 2 * p3.XYZ()) / 6;
  const gp_XYZ b2 = (2 * p0.XYZ() - 9 * p1.XYZ() + 18 * p2.XYZ() - 5 * p3.XYZ()) / 6;
  NCollection_Array1<gp_Pnt> poles(1, 4);
  poles(1) = p0; poles(2) = gp_Pnt(b1); poles(3) = gp_Pnt(b2); poles(4) = p3;
  return new Geom_BezierCurve(poles);
}

void patch(double a) {
  const TopoDS_Edge e0 = BRepBuilderAPI_MakeEdge(side(a, true, 0)).Edge();
  const TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(side(a, false, 1)).Edge();
  Handle(Geom_BezierCurve) c2 = side(a, true, 1); c2->Reverse();
  const TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(c2).Edge();
  Handle(Geom_BezierCurve) c3 = side(a, false, 0); c3->Reverse();
  const TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(c3).Edge();
  char input[64];
  std::snprintf(input, sizeof input, "cubic quad, a=%.2f", a);

  const auto quality = [&](const TopoDS_Shape& made, double& interior, double& boundary) {
    interior = boundary = 0;
    TopoDS_Face face;
    for (TopExp_Explorer it(made, TopAbs_FACE); it.More(); it.Next()) face = TopoDS::Face(it.Current());
    if (face.IsNull()) { interior = boundary = -1; return; }
    BRepAdaptor_Surface surface(face);
    for (int i = 0; i <= 32; ++i) for (int j = 0; j <= 32; ++j) {
      const gp_Pnt expected = surfacePoint(i / 32., j / 32., a);
      GeomAPI_ProjectPointOnSurf project(expected, BRep_Tool::Surface(face));
      interior = std::max(interior, project.NbPoints() ? project.LowerDistance() : 1.0);
    }
    for (const TopoDS_Edge& edge : {e0, e1, e2, e3}) {
      BRepAdaptor_Curve curve(edge);
      for (int k = 0; k <= 128; ++k) {
        const gp_Pnt p = curve.Value(curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * k / 128.);
        GeomAPI_ProjectPointOnSurf project(p, BRep_Tool::Surface(face));
        boundary = std::max(boundary, project.NbPoints() ? project.LowerDistance() : 1.0);
      }
    }
  };
  TopoDS_Shape stockResult;
  const Stats stock = measure([&] {
    BRepFill_Filling filling(3, 15, 2, false, 1e-5, 1e-4, 1e-2, 0.1, 8, 9);
    for (const TopoDS_Edge& edge : {e0, e1, e2, e3}) filling.Add(edge, GeomAbs_C0);
    filling.Build();
    if (!filling.IsDone()) return false;
    stockResult = filling.Face();
    return true;
  });
  TopoDS_Shape extResult;
  const std::vector<uint32_t> edges{BRepExt::store(session(), e0), BRepExt::store(session(), e1), BRepExt::store(session(), e2), BRepExt::store(session(), e3)};
  const Stats ext = measure([&] {
    const auto result = BRepExt::patch(session(), edges, std::vector<uint32_t>(4, 0), std::vector<int>(4, 0), {}, {}, 1e-6, 1, 0.1, 0);
    extResult = BRepExt::get(session(), result.evolution.resultId);
    return !extResult.IsNull();
  });
  double interior, boundary;
  char q[128];
  quality(stockResult, interior, boundary);
  std::snprintf(q, sizeof q, "off surface: interior %.3e, sides %.3e", interior, boundary);
  row("patch", input, "stock", stock, q);
  quality(extResult, interior, boundary);
  std::snprintf(q, sizeof q, "off surface: interior %.3e, sides %.3e", interior, boundary);
  row("patch", input, "extension", ext, q);
}

// ---- hollow ---------------------------------------------------------------------------------
// A box a x b x c with one edge rounded by r, hollowed to a wall t with the outer skin kept.
// The cavity is the inward offset by t: the box (a - 2t)(b - 2t)(c - 2t), its rounded edge having
// radius r - t where t < r (volume less (1 - pi/4)(r - t)^2 (c - 2t)) and being sharp for t >= r.
// The wall holds V_outer - V_cavity with V_outer = abc - (1 - pi/4) r^2 c.
void hollow(double r, double t) {
  const double a = 60, b = 40, c = 30;
  TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  TopoDS_Edge edge;
  for (TopExp_Explorer it(box, TopAbs_EDGE); it.More(); it.Next()) {
    const TopoDS_Edge e = TopoDS::Edge(it.Current());
    BRepAdaptor_Curve curve(e);
    const gp_Pnt p0 = curve.Value(curve.FirstParameter()), p1 = curve.Value(curve.LastParameter());
    if (std::abs(p0.X() - a) < 1e-9 && std::abs(p1.X() - a) < 1e-9 && std::abs(p0.Y() - b) < 1e-9 && std::abs(p1.Y() - b) < 1e-9) edge = e;
  }
  BRepFilletAPI_MakeFillet maker(box);
  maker.Add(r, edge);
  maker.Build();
  TopoDS_Shape body = maker.Shape();
  for (TopExp_Explorer it(maker.Shape(), TopAbs_SOLID); it.More(); it.Next()) body = it.Current();
  const double outer = a * b * c - (1 - M_PI / 4) * r * r * c;
  const double cavity = (a - 2 * t) * (b - 2 * t) * (c - 2 * t) - (t < r ? (1 - M_PI / 4) * (r - t) * (r - t) * (c - 2 * t) : 0);
  const double expected = outer - cavity;
  char input[64];
  std::snprintf(input, sizeof input, "box, round r=%.1f, wall t=%.1f", r, t);
  double errorStock = -1, errorExt = -1;
  bool validStock = false, validExt = false;
  const Stats stock = measure([&] {
    BRepOffsetAPI_MakeThickSolid thick;
    thick.MakeThickSolidByJoin(body, NCollection_List<TopoDS_Shape>(), -t, 1e-4);
    if (!thick.IsDone()) return false;
    validStock = valid(thick.Shape());
    errorStock = std::abs(volume(thick.Shape()) - expected);
    return true;
  });
  const uint32_t bodyId = BRepExt::store(session(), body);
  const Stats ext = measure([&] {
    const auto result = BRepExt::hollow(session(), bodyId, t, 1e-4, {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    validExt = valid(made);
    errorExt = std::abs(volume(made) - expected);
    return true;
  });
  char q[96];
  std::snprintf(q, sizeof q, "%s, |dV| %.3e", validStock ? "valid" : "INVALID", errorStock);
  row("hollow", input, "stock", stock, q);
  std::snprintf(q, sizeof q, "%s, |dV| %.3e", validExt ? "valid" : "INVALID", errorExt);
  row("hollow", input, "extension", ext, q);
}

// ---- offset ---------------------------------------------------------------------------------
void offset(double d) {
  const double a = 60, b = 40;
  const gp_Pnt p[4] = {gp_Pnt(0, 0, 0), gp_Pnt(a, 0, 0), gp_Pnt(a, b, 0), gp_Pnt(0, b, 0)};
  BRepBuilderAPI_MakeWire maker;
  for (int i = 0; i < 4; ++i) maker.Add(BRepBuilderAPI_MakeEdge(p[i], p[(i + 1) % 4]).Edge());
  const TopoDS_Wire wire = maker.Wire();
  const double expected = 2 * (a + b) + 2 * M_PI * d;
  char input[64];
  std::snprintf(input, sizeof input, "rectangle, outward d=%.1f", d);
  double errorStock = -1, errorExt = -1;
  const Stats stock = measure([&] {
    BRepOffsetAPI_MakeOffset offset(wire, GeomAbs_Arc);
    offset.Perform(d);
    if (!offset.IsDone()) return false;
    errorStock = std::abs(length(offset.Shape()) - expected);
    return true;
  });
  const uint32_t id = BRepExt::store(session(), wire);
  const Stats ext = measure([&] {
    const auto result = BRepExt::offsetCurve(session(), id, d, 0);
    errorExt = std::abs(length(BRepExt::get(session(), result.resultId)) - expected);
    return true;
  });
  char q[96];
  std::snprintf(q, sizeof q, "|dL| %.3e", errorStock);
  row("offset", input, "stock", stock, q);
  std::snprintf(q, sizeof q, "|dL| %.3e", errorExt);
  row("offset", input, "extension", ext, q);
}

// ---- rounds beyond the easy case ------------------------------------------------------------
// Several edges at once, rounds that run into earlier rounds, and corners where blends meet. A body
// is a box a x b x c (40 x 30 x 20). `all` has the closed form of Steiner's formula for a box
// rounded by r everywhere: with x = a - 2r etc., V = xyz + 2r(xy + yz + zx) + pi r^2 (x + y + z) +
// (4/3) pi r^3 (the Minkowski sum of the inner box and a ball). The other rows have no closed form and
// report validity and volume only.
std::vector<TopoDS_Edge> edgesWhere(const TopoDS_Shape& shape, const std::function<bool(const Bnd_Box&, double*)>& pick) {
  std::vector<TopoDS_Edge> out;
  for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next()) {
    Bnd_Box box;
    BRepBndLib::Add(it.Current(), box);
    double b[6];
    box.Get(b[0], b[1], b[2], b[3], b[4], b[5]);
    if (pick(box, b)) out.push_back(TopoDS::Edge(it.Current()));
  }
  return out;
}
bool near(double a, double b) { return std::abs(a - b) < 1e-6; }
using Pick = std::function<bool(const Bnd_Box&, double*)>;
const Pick everyEdge = [](const Bnd_Box&, double*) { return true; };
const Pick topEdges = [](const Bnd_Box&, double* b) { return near(b[2], 20) && near(b[5], 20); };
const Pick uprights = [](const Bnd_Box&, double* b) { return near(b[0], b[3]) && near(b[1], b[4]) && b[5] - b[2] > 5; };
const Pick cornerEdges = [](const Bnd_Box&, double* b) {
  const bool atX = near(b[0], 0) && near(b[3], 0), atY = near(b[1], 0) && near(b[4], 0);
  return b[0] <= 1e-6 && b[1] <= 1e-6 && b[5] >= 20 - 1e-6 && (atX || atY || (near(b[0], b[3]) && near(b[1], b[4])));
};

struct Step { Pick pick; double radius; };
void rounds(const char* name, const std::vector<Step>& steps, double expected) {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(40, 30, 20).Shape();
  double volumeStock = 0, volumeExt = 0;
  bool validStock = false, validExt = false;
  const Stats stock = measure([&] {
    TopoDS_Shape body = box;
    for (const Step& step : steps) {
      BRepFilletAPI_MakeFillet maker(body);
      for (const TopoDS_Edge& edge : edgesWhere(body, step.pick)) maker.Add(step.radius, edge);
      maker.Build();
      if (!maker.IsDone()) return false;
      body = maker.Shape();
    }
    validStock = valid(body);
    volumeStock = volume(body);
    return validStock;
  });
  const Stats checked = measure([&] {
    TopoDS_Shape body = box;
    for (const Step& step : steps) {
      BRepFilletAPI_MakeFillet maker(body);
      for (const TopoDS_Edge& edge : edgesWhere(body, step.pick)) maker.Add(step.radius, edge);
      maker.Build();
      if (!maker.IsDone()) return false;
      body = maker.Shape();
    }
    return validTimed(body);
  });
  const Stats ext = measure([&] {
    TopoDS_Shape body = box;
    for (const Step& step : steps) {
      BRepExt::blend2::Request request;
      request.operation = "round";
      request.body = body;
      for (const TopoDS_Edge& edge : edgesWhere(body, step.pick)) request.edges.push_back(edge);
      request.first = step.radius;
      request.law = BRepExt::blend2::Law::Conic;
      request.rho = 0.5;
      body = BRepExt::blend2::blend(request).shape;
    }
    validExt = valid(body);
    volumeExt = volume(body);
    return validExt;
  });
  char q[200];
  const auto text = [&](bool ok, double v) {
    if (expected > 0) std::snprintf(q, sizeof q, "%s, V %.6f, |dV| %.3e", ok ? "valid" : "INVALID", v, std::abs(v - expected));
    else std::snprintf(q, sizeof q, "%s, V %.6f", ok ? "valid" : "INVALID", v);
  };
  text(validStock, volumeStock);
  row("round", name, "stock", stock, q);
  row("round", name, "stock+check", checked, q);
  text(validExt, volumeExt);
  row("round", name, "extension", ext, q);
}
// Rounds along curved junctions: two crossing cylinders and a sphere standing on a box. The edges are
// space curves where two curved surfaces meet, the case where a blend must follow a surface-surface
// intersection and change its contact curves as it goes. No closed form: validity and volume only,
// and the volumes of the two implementations are compared.
void junction(const char* name, const TopoDS_Shape& body, const Pick& pick, double radius) {
  double volumeStock = 0, volumeExt = 0;
  bool validStock = false, validExt = false;
  const Stats stock = measure([&] {
    BRepFilletAPI_MakeFillet maker(body);
    for (const TopoDS_Edge& edge : edgesWhere(body, pick)) maker.Add(radius, edge);
    maker.Build();
    if (!maker.IsDone()) return false;
    validStock = valid(maker.Shape());
    volumeStock = volume(maker.Shape());
    return validStock;
  });
  const Stats checked = measure([&] {
    BRepFilletAPI_MakeFillet maker(body);
    for (const TopoDS_Edge& edge : edgesWhere(body, pick)) maker.Add(radius, edge);
    maker.Build();
    return maker.IsDone() && validTimed(maker.Shape());
  });
  const Stats ext = measure([&] {
    BRepExt::blend2::Request request;
    request.operation = "round";
    request.body = body;
    for (const TopoDS_Edge& edge : edgesWhere(body, pick)) request.edges.push_back(edge);
    request.first = radius;
    request.law = BRepExt::blend2::Law::Conic;
    request.rho = 0.5;
    const TopoDS_Shape made = BRepExt::blend2::blend(request).shape;
    validExt = valid(made);
    volumeExt = volume(made);
    return validExt;
  });
  char q[160];
  std::snprintf(q, sizeof q, "%s, V %.6f", validStock ? "valid" : "INVALID", volumeStock);
  row("junction", name, "stock", stock, q);
  row("junction", name, "stock+check", checked, q);
  std::snprintf(q, sizeof q, "%s, V %.6f", validExt ? "valid" : "INVALID", volumeExt);
  row("junction", name, "extension", ext, q);
}

double steiner(double r) {
  const double x = 40 - 2 * r, y = 30 - 2 * r, z = 20 - 2 * r;
  return x * y * z + 2 * r * (x * y + y * z + z * x) + M_PI * r * r * (x + y + z) + 4. / 3 * M_PI * r * r * r;
}

} // namespace

int main(int argc, char** argv) {
  // With an argument, only that group is run (for debugging one group: fillet, patch, rounds, junction, hollow or offset).
  const std::string only = argc > 1 ? argv[1] : "";
  const auto want = [&](const char* group) { return only.empty() || only == group; };
  std::printf("operation | input | who | ok/runs | CPU median [min, max] ms | quality\n");
  if (want("fillet")) for (double r : {0.2, 0.5}) fillet(r);
  if (want("patch")) for (double a : {0.5, 2.0}) patch(a);
  if (want("rounds")) {
    rounds("box, all 12 edges r=6", {{everyEdge, 6}}, steiner(6));
    rounds("box, all 12 edges r=9 (2r near height)", {{everyEdge, 9}}, steiner(9));
    rounds("box, three edges at a corner r=7", {{cornerEdges, 7}}, 0);
    rounds("box, top r=8 then uprights r=3", {{topEdges, 8}, {uprights, 3}}, 0);
  }
  if (want("junction")) {
    const TopoDS_Shape along = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -20), gp_Dir(0, 0, 1)), 10, 40).Shape();
    const TopoDS_Shape across = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-20, 0, 0), gp_Dir(1, 0, 0)), 6, 40).Shape();
    const TopoDS_Shape crossed = BRepAlgoAPI_Fuse(along, across).Shape();
    const Pick inside = [](const Bnd_Box&, double* b) { return b[3] - b[0] < 39 && b[5] - b[2] < 39 && b[3] - b[0] > 1; };
    for (const double r : {1.0, 2.5}) {
      char label[64];
      std::snprintf(label, sizeof label, "crossed cylinders R10 and R6, r=%.1f", r);
      junction(label, crossed, inside, r);
    }
    const TopoDS_Shape standing = BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(gp_Pnt(-25, -25, -10), 50, 50, 10).Shape(),
                                                   BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 14).Shape()).Shape();
    const Pick ring = [](const Bnd_Box&, double* b) { return near(b[2], b[5]) && near(b[2], 0) && b[3] - b[0] < 29 && b[3] - b[0] > 20; };
    for (const double r : {2.0, 6.0}) {
      char label[64];
      std::snprintf(label, sizeof label, "sphere R14 on a box, r=%.1f", r);
      junction(label, standing, ring, r);
    }
  }
  if (want("hollow")) {
    hollow(3, 2);
    hollow(3, 4);
  }
  if (want("offset")) for (double d : {2.0, 10.0}) offset(d);
  return 0;
}
