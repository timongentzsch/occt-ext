// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <BRepBlendExt_Internal.hxx>
#include <GeomAPI_IntSS.hxx>
#include <Geom_OffsetSurface.hxx>

namespace BRepExt::blending::detail
{

bool intersectSpine(Rail& rail)
{
  // Closed, constant-distance circuits have no stop or size-law correspondence to transfer.
  // Open chains need endpoint correspondence before they can use an independent carrier.
  if (!rail.closed || !rail.sizes.places.empty() || rail.angle != 0
      || (rail.law != Law::Conic && rail.law != Law::G2)
      || rail.from != rail.first || rail.to != rail.last)
    return false;
  for (const Side& side : rail.s)
    if (side.cliffCurve || !side.trim.IsNull())
      return false;

  const double middle = (rail.first + rail.last) / 2;
  Station seed;
  seed.t = middle;
  for (int i = 0; i < 2; ++i)
  {
    double first, last;
    const auto pc = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, first, last);
    if (pc.IsNull())
      return false;
    const gp_Pnt2d uv = pc->Value(middle);
    seed.uv[i][0] = uv.X();
    seed.uv[i][1] = uv.Y();
  }
  // Grow a seed on the original branch, without asking its station map to traverse a circuit.
  for (int i = 1; i <= 8; ++i)
    if (!contactAt(rail, middle, seed, i / 8.0))
      return false;

  Handle(Geom_Surface) offsets[2];
  for (int i = 0; i < 2; ++i)
  {
    const Handle(Geom_OffsetSurface) offset = new Geom_OffsetSurface(rail.s[i].surface,
      rail.sigma * rail.s[i].distance * rail.s[i].outward);
    // Retain exact canonical offsets where the geometry library supplies them; otherwise
    // the same intersector receives the general offset evaluator.
    offsets[i] = offset->Surface();
    if (offsets[i].IsNull())
      offsets[i] = offset;
  }
  // Use the general intersector, including its analytic and marching cases. No support-family
  // formula is used here. The resulting curve is only a carrier: ballAt corrects on both
  // original offsets at every station, so an approximate SSI curve is not a contact curve.
  GeomAPI_IntSS intersection(offsets[0], offsets[1], rail.resolved(limit::solved));
  if (!intersection.IsDone())
    return false;
  Handle(Geom_Curve) chosen;
  double at = 0;
  for (int i = 1; i <= intersection.NbLines(); ++i)
  {
    const auto curve = intersection.Line(i);
    if (!curve->IsClosed())
      continue;
    GeomAPI_ProjectPointOnCurve projection(seed.centre, curve,
                                           curve->FirstParameter(), curve->LastParameter());
    if (projection.NbPoints() == 0
        || projection.LowerDistance() > rail.resolved(limit::together))
      continue;
    if (!chosen.IsNull())
      return false; // Ambiguous branch identification is not a continuation certificate.
    chosen = curve;
    at = projection.LowerDistanceParameter();
  }
  if (chosen.IsNull())
    return false;
  gp_Pnt point;
  gp_Vec tangent;
  chosen->D1(at, point, tangent);
  const double orientation = tangent.Dot(seed.along);
  if (!(tangent.Magnitude() > 0) || orientation == 0)
    return false;
  gp_Pnt start, end;
  gp_Vec ds, de;
  chosen->D1(chosen->FirstParameter(), start, ds);
  chosen->D1(chosen->LastParameter(), end, de);
  if (start.Distance(end) > rail.resolved(limit::solved)
      || !(ds.Magnitude() > 0 && de.Magnitude() > 0)
      || ds.Normalized().Dot(de.Normalized()) < 1 - tol::angle)
    return false;
  const double span = (chosen->LastParameter() - chosen->FirstParameter())
                      * (orientation < 0 ? -1 : 1);
  rail.spineParameterFirst = rail.first;
  rail.spineParameterLast  = rail.last;
  rail.spine               = chosen;
  rail.spineFirst          = at - span / 2;
  rail.spineLast           = at + span / 2;
  return true;
}

//! Analytic construction for a cylindrical cavity meeting a transverse plane.
//! The circular envelope crosses when r >= (R + r)|n_plane.axis|. The window
//! h=2 asin(r/(R + r)) encloses the complete ellipse double-point domain.
//! This is a construction domain, not an approximation allowance.
bool prepareConicWindow(Rail& rail)
{
  if (rail.spine.IsNull() || rail.law != Law::Conic || rail.rho != 0.5 || rail.sigma != -1)
    return false;
  for (int c = 0; c < 2; ++c)
  {
    GeomAdaptor_Surface cylinder(rail.s[c].surface), plane(rail.s[1 - c].surface);
    if (cylinder.GetType() != GeomAbs_Cylinder || plane.GetType() != GeomAbs_Plane
        || rail.s[c].outward != -1 || rail.s[c].distance != rail.s[1 - c].distance)
      continue;
    const gp_Vec axis(cylinder.Cylinder().Axis().Direction());
    gp_Vec normal(plane.Plane().Axis().Direction());
    normal *= rail.s[1 - c].outward;
    const double cosine = normal.Dot(axis), r = rail.s[c].distance;
    const double R = cylinder.Cylinder().Radius();
    if (!(r > 0 && r < R && std::abs(cosine) > tol::angle
          && r >= (R + r) * std::abs(cosine)))
      return false;
    gp_Vec acute = normal - axis * cosine;
    if (!(acute.Magnitude() > tol::angle))
      return false;
    acute.Normalize();
    rail.conicWindow = Rail::ConicWindow{c,
                                        axis,
                                        acute,
                                        normal,
                                        plane.Plane().Location(),
                                        2 * std::asin(r / (R + r))};
    return true;
  }
  return false;
}

double conicPhase(const Rail& rail, const Station& station)
{
  const auto& window = *rail.conicWindow;
  const gp_Vec radial = station.normal[window.cylinder] * rail.s[window.cylinder].outward;
  return std::atan2(window.axis.Dot(window.acute.Crossed(radial)), window.acute.Dot(radial));
}

} // namespace BRepExt::blending::detail
