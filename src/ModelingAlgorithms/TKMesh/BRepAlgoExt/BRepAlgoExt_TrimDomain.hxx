// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepTopAdaptor_FClass2d.hxx>
#include <Bnd_Box2d.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <memory>
#include <utility>
#include <vector>

namespace BRepExt
{
// Endpoints and stationary parameters of one rational Bezier coordinate, using the same
// Bernstein root isolation as the trim's monotone pcurve pieces.
std::vector<double> coordinateExtrema(const Handle(Geom_BezierCurve)&, int axis);

// The trimmed part of a face's parameter plane, from its pcurves: what trim-only area,
// centroid and bounds integrate or refine over, never the rest of the support surface.
class TrimDomain
{
  TopoDS_Face             face;
  BRepTopAdaptor_FClass2d classifier;

  // A run of one pcurve along which v is monotone, across the junctions of its Bezier arcs:
  // a row of constant v crosses it once at most.
  struct Segment
  {
    Handle(Geom2d_BezierCurve) arc;
    double                     t0, t1, y0, y1;
  };

  struct Piece
  {
    std::vector<Segment> segments;
    double               y0, y1;
  };

  std::vector<Piece> pieces;

  // Pole boxes of the boundary arcs, halved on demand to tell cells it crosses.
  struct Boundary
  {
    Handle(Geom2d_BezierCurve)        curve;
    Bnd_Box2d                         box;
    mutable std::unique_ptr<Boundary> left, right;
    explicit Boundary(const Handle(Geom2d_BezierCurve)&);
    bool intersects(double, double, double, double, int = 0) const;
  };

  std::vector<std::unique_ptr<Boundary>> boundaries;
  double                                 crossing(const Piece&, double v) const;

public:
  // A band of rows between consecutive breaks: the same pieces cross every row of it, in
  // the same order, with the same gaps between them inside the trim.
  struct Band
  {
    double            lo, hi;
    std::vector<int>  order;
    std::vector<char> inside;
  };

  double u0, u1, v0, v1;
  explicit TrimDomain(const TopoDS_Face&);
  bool              boundaryIntersects(double a, double b, double c, double d) const;
  std::vector<Band> bands(std::vector<double> breaks) const;
  void              spans(const Band&, double v, std::vector<std::pair<double, double>>& out) const;
  // Whether a point is inside, from the bands' crossings of its row.
  bool inside(const std::vector<Band>&, double u, double v) const;
};

struct TrimProperties
{
  double area = 0;
  gp_Pnt center;
  double volume    = 0;
  bool   converged = false;
};

// Area and centroid integrated over the trim only; optionally the forward face's
// signed volume contribution about an origin. Throws where the trim cannot be read.
TrimProperties trimmedProperties(const TopoDS_Face&, const gp_Pnt* volumeOrigin = nullptr);
} // namespace BRepExt
