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
    bool ok = false;
    try { ok = work(); }
    catch (const Standard_Failure& e) { ok = false; s.refusal = e.GetMessageString() ? e.GetMessageString() : "Standard_Failure"; }
    catch (const std::exception& e) { ok = false; s.refusal = e.what(); }
    const double ms = 1000. * (std::clock() - begin) / CLOCKS_PER_SEC;
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

double volume(const TopoDS_Shape& shape) {
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  return props.Mass();
}
double length(const TopoDS_Shape& shape) {
  GProp_GProps props;
  BRepGProp::LinearProperties(shape, props);
  return props.Mass();
}
bool valid(const TopoDS_Shape& shape) { return !shape.IsNull() && BRepCheck_Analyzer(shape).IsValid(); }

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

} // namespace

int main() {
  std::printf("operation | input | who | ok/runs | CPU median [min, max] ms | quality\n");
  for (double r : {0.2, 0.5}) fillet(r);
  for (double a : {0.5, 2.0}) patch(a);
  hollow(3, 2);
  hollow(3, 4);
  for (double d : {2.0, 10.0}) offset(d);
  return 0;
}
