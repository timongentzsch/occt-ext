// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Spiral: a helix with a spiral pitch. At t turns of [0, turns] the wire is at angle
// 2πt about the axis, `radius` + `pitch`·t from it and `height`·t/turns along it: the radius
// grows by the pitch each turn while the rise stays linear.

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>

#include <AdvApprox_ApproxAFunction.hxx>
#include <AdvApprox_EvaluatorFunction.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Wire.hxx>

#include <cmath>

namespace BRepExt
{
namespace
{

// The spiral as a function of turns t. Every derivative is exact, so the approximation is held
// to the tolerance by the curve itself, not by samples of it.
struct SpiralLaw final : AdvApprox_EvaluatorFunction
{
  gp_XYZ origin, x, y, z;
  double radius = 0, pitch = 0, rise = 0;

  void Evaluate(int* dimension,
                double /*startEnd*/[2],
                double* parameter,
                int*    order,
                double* result,
                int*    error) override
  {
    if (*dimension != 3)
    {
      *error = 1;
      return;
    }
    const double t = *parameter, angle = 2 * M_PI * t;
    // The m-th derivative of the unit radial cos(2πt)·x + sin(2πt)·y.
    auto radial = [&](int m) {
      const double scale = std::pow(2 * M_PI, m), phase = angle + m * M_PI / 2;
      return (x * std::cos(phase) + y * std::sin(phase)) * scale;
    };
    // Leibniz on (radius + pitch·t)·radial(t): pitch is the only derivative of the radius.
    gp_XYZ value = radial(*order) * (radius + pitch * t);
    if (*order >= 1)
      value += radial(*order - 1) * (*order * pitch);
    if (*order == 0)
      value += origin + z * (rise * t);
    if (*order == 1)
      value += z * rise;
    result[0] = value.X();
    result[1] = value.Y();
    result[2] = value.Z();
    *error    = 0;
  }
};

} // namespace

EvolutionData spiral(OcctKernel& kernel,
                     Vec3        origin,
                     Vec3        axis,
                     Vec3        start,
                     double      radius,
                     double      height,
                     double      turns,
                     double      pitch,
                     bool        leftHanded,
                     double      tolerance)
{
  const char* op = "spiral";
  if (!finite(origin) || !finite(axis) || !finite(start) || !std::isfinite(radius)
      || !std::isfinite(height) || !std::isfinite(turns) || !std::isfinite(pitch)
      || !std::isfinite(tolerance))
    fail(op, invalid_input, "inputs must be finite");
  if (!(turns > 0))
    fail(op, invalid_input, "turns must be positive");
  if (radius < 0 || height < 0)
    fail(op, invalid_input, "radius and height must not be negative");
  if (!(tolerance > 0))
    fail(op, invalid_input, "tolerance must be positive");
  if (radius + pitch * turns < 0)
    fail(op, invalid_input, "the spiral would pass through its axis");
  if (height == 0 && pitch == 0)
    fail(op, invalid_input, "a spiral without length or pitch retraces a circle");
  if (radius == 0 && pitch == 0)
    fail(op, invalid_input, "a spiral without radius or pitch is its axis");
  const gp_XYZ z(axis.x, axis.y, axis.z), along(start.x, start.y, start.z);
  if (z.Modulus() == 0)
    fail(op, invalid_input, "the axis must have a direction");
  const gp_XYZ unitZ = z / z.Modulus(), radial = along - unitZ * along.Dot(unitZ);
  if (radial.Modulus() <= BRepExt::tolerance::angle * along.Modulus())
    fail(op, invalid_input, "the start direction must not lie along the axis");
  SpiralLaw law;
  law.origin = gp_XYZ(origin.x, origin.y, origin.z);
  law.z      = unitZ;
  law.x      = radial / radial.Modulus();
  law.y      = unitZ.Crossed(law.x) * (leftHanded ? -1.0 : 1.0);
  law.radius = radius;
  law.pitch  = pitch;
  law.rise   = height / turns;
  try
  {
    // C2 needs degree 5 at least; seven leaves room for few spans. A span of a 64th of a
    // turn is far finer than either needs, so it bounds the search.
    const int                                 maxSpans = static_cast<int>(std::ceil(turns)) * 64;
    const Handle(NCollection_HArray1<double>) none,
      threeD = new NCollection_HArray1<double>(1, 1, tolerance);
    AdvApprox_ApproxAFunction
      fit(0, 0, 1, none, none, threeD, 0, turns, GeomAbs_C2, 7, maxSpans, law);
    if (!fit.IsDone() || !fit.HasResult())
      fail(op, no_solution, "the spiral could not be approximated");
    if (fit.MaxError(3, 1) > tolerance)
      fail(op, no_solution, "the spiral's approximation exceeds the tolerance");
    NCollection_Array1<gp_Pnt> poles(1, fit.NbPoles());
    fit.Poles(1, poles);
    Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles,
                                                            fit.Knots()->Array1(),
                                                            fit.Multiplicities()->Array1(),
                                                            fit.Degree());
    BRepBuilderAPI_MakeWire   wire(BRepBuilderAPI_MakeEdge(curve).Edge());
    if (!wire.IsDone())
      fail(op, no_solution, "the spiral's wire could not be made");
    EvolutionData out;
    out.resultId = store(kernel, wire.Wire());
    return out;
  }
  catch (const Standard_Failure& e)
  {
    fail(op, no_solution, e.what());
  }
}

} // namespace BRepExt
