// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepBlendExt_Request.hxx>
#include <BRepToolsExt_Session.hxx>

// The operations this package adds to the kernel. Each takes and returns shape handles (see
// BRepToolsExt_Session.hxx), reports face history like upstream's *WithHistory methods (see
// BRepToolsExt_History.hxx) and fails through BRepExt::fail (see GeomAlgoExt_Failure.hxx).
// Lengths are in model units. `hashes` and `bound` name the input faces whose history is
// tracked (BRepToolsExt_History.hxx). The numeric values of the enumerations are fixed.

namespace BRepExt
{
struct Failure;

//! How the faces of an offset are joined where they part; the values are GeomAbs_JoinType's.
enum class OffsetJoin
{
  Arc          = 0,
  Tangent      = 1,
  Intersection = 2
};

//! The boundary of a rectangular sheet in its parameters.
enum class SheetBoundary
{
  UFirst = 0,
  ULast  = 1,
  VFirst = 2,
  VLast  = 3
};

//! The shape a curve is extended by at an end.
enum class CurveExtension
{
  Natural    = 0, //!< the curve's own continuation; conics stay conics
  Linear     = 1, //!< a line along the end tangent
  Soft       = 2, //!< G2 at the join, curvature falling to zero at the far end
  Reflective = 3, //!< the end segment mirrored in the plane square to the curve at its end
  Arc        = 4  //!< the circle of the end curvature
};

//! What `rebuildCurve` holds to.
enum class CurveRebuild
{
  Tolerance = 0, //!< approximated within the tolerance
  Count     = 1  //!< a fixed degree and count of poles
};

//! A shape-preserving curve edit.
enum class CurveRefinement
{
  RaiseDegree = 0, //!< each spline's degree raised by the amount
  Subdivide   = 1, //!< each spline's poles doubled by knot insertion
  Simplify    = 2  //!< redundant edges merged within the tolerance
};

//! How the gap at a corner of a curve offset is closed.
enum class OffsetGap
{
  Round   = 0, //!< by an arc about the corner
  Linear  = 1, //!< by the chord of that arc
  Natural = 2  //!< by carrying the two neighbours on to where they cross
};

//! A recognised primitive.
enum class PrimitiveKind
{
  None     = 0,
  Box      = 1,
  Cylinder = 2,
  Sphere   = 3
};

//! For one flat face pulled out to a wall that faces the same way: what the face does there.
enum class FaceGrow
{
  None   = 0, //!< goes on past the wall
  Fixed  = 1, //!< stops flush with the wall
  Moving = 2  //!< goes on with the wall
};

//! How two guides place a swept profile.
enum class GuideMethod
{
  Point = 0, //!< turned and scaled so its point on the second guide stays on it
  Chord = 1, //!< only turned
  Curve = 2  //!< turned so that it passes through the second guide
};

//! How curves are imprinted on a body.
enum class ImprintMethod
{
  Direction = 0, //!< along a direction
  Normal    = 1  //!< along the normals of the faces
};

//! What becomes of the loose ends of imprinted edges.
enum class ImprintCompletion
{
  None   = 0, //!< left as they are
  Extend = 1  //!< carried on straight in the face's parameters until they leave it
};

//! The corners of a swept section; the values are BRepBuilderAPI_TransitionMode's.
enum class SweepCorner
{
  Transformed = 0,
  Right       = 1,
  Round       = 2
};

//! Clear the triangulation of the shape behind `shape`. Defined by the embedding program.
Standard_EXPORT void clearTriangulation(OcctKernel& kernel, uint32_t shape);

//! An axis-aligned box from the origin with sides `dx`, `dy`, `dz`. Fails where a side is at
//! most the kernel resolution.
Standard_EXPORT EvolutionData makeBox(OcctKernel& kernel, double dx, double dy, double dz);

//! The union, difference or intersection (`operation`) of `sources` with `tools`, verified as
//! booleanCells; faces within `angle` (radians, 0: none) are merged. Fails on another selection.
Standard_EXPORT EvolutionData booleanChecked(OcctKernel&                  kernel,
                                             Select                       operation,
                                             const std::vector<uint32_t>& sources,
                                             const std::vector<uint32_t>& tools,
                                             BooleanGlue                  glue,
                                             double                       fuzzy,
                                             double                       angle,
                                             const std::vector<int>&      hashes,
                                             int                          bound);

//! The cells of one General Fuse that `select` keeps, Union to Region (BRepAlgoExt_Cells.hxx).
//! Fails where no attempt's partition is verified, or on a later selection.
Standard_EXPORT EvolutionData booleanCells(OcctKernel&                  kernel,
                                           Select                       select,
                                           const std::vector<uint32_t>& sources,
                                           const std::vector<uint32_t>& tools,
                                           double                       fuzzy,
                                           const std::vector<int>&      hashes,
                                           int                          bound);

//! What each attempt of booleanCells did for that selection, whether or not one was accepted.
//! Fails only where the request is refused before any attempt.
Standard_EXPORT std::vector<BooleanAttempt> booleanReport(OcctKernel&                  kernel,
                                                          Select                       select,
                                                          const std::vector<uint32_t>& sources,
                                                          const std::vector<uint32_t>& tools,
                                                          double                       fuzzy);

//! A solid hollowed to a wall of `thickness` with `openFaces` removed; with `local`, only the wall
//! round the openings. Fails where the wall is within the clearance allowance of `tolerance`.
Standard_EXPORT EvolutionData shell(OcctKernel&                  kernel,
                                    uint32_t                     body,
                                    const std::vector<uint32_t>& openFaces,
                                    double                       thickness,
                                    double                       tolerance,
                                    OffsetJoin                   join,
                                    const std::vector<int>&      hashes,
                                    int                          bound,
                                    bool                         local);

//! A B-spline face's degrees and span counts: U degree, V degree, U spans, V spans. Fails on
//! another surface.
Standard_EXPORT std::vector<int> surfaceSplineInfo(OcctKernel& kernel, uint32_t face);
struct Vec3;
//! A closed periodic cubic spline edge through `coordinates` (x, y, z in turn), with chord-length
//! knots. Fails for fewer than three points or repeated successive points.
Standard_EXPORT EvolutionData interpolateClosed(OcctKernel&                kernel,
                                                const std::vector<double>& coordinates);
//! A face closing round through three or more open single-edge `sections`. Fails where the
//! periodic interpolation misses `tolerance`.
Standard_EXPORT EvolutionData skinPeriodic(OcctKernel&                  kernel,
                                           const std::vector<uint32_t>& sections,
                                           double                       tolerance);
//! The shapes of IGES file content. Defined by the embedding program.
Standard_EXPORT EvolutionData importIges(OcctKernel& kernel, const std::string& text);
//! The shape behind `shape` as IGES file content. Defined by the embedding program.
Standard_EXPORT std::string exportIges(OcctKernel& kernel, uint32_t shape);
//! A lone face whose support is refitted as a B-spline of the degrees and `spans`, its trims and
//! holes kept. Fails where the trim boundary misses `tolerance`.
Standard_EXPORT EvolutionData rebuildFace(OcctKernel& kernel,
                                          uint32_t    face,
                                          int         uDegree,
                                          int         vDegree,
                                          int         spans,
                                          double      tolerance);
//! Rebuild one face of a body in place: its support is refitted on the face's UV box and its
//! edges stay shared with the neighbours, within the fit tolerance.
Standard_EXPORT EvolutionData rebuildBodyFace(OcctKernel&             kernel,
                                              uint32_t                body,
                                              uint32_t                face,
                                              int                     uDegree,
                                              int                     vDegree,
                                              int                     spans,
                                              double                  tolerance,
                                              const std::vector<int>& hashes,
                                              int                     bound);
//! A B-spline or Bezier face support as exact world-space B-spline data: degrees,
//! knot and pole counts, knots with
//! multiplicities in U then V, then each pole as x, y, z, weight, U-major.
Standard_EXPORT std::vector<double> faceSpline(OcctKernel& kernel, uint32_t face);
//! A face on the whole of the spline that `faceSpline` describes.
Standard_EXPORT EvolutionData splineFace(OcctKernel& kernel, const std::vector<double>& spline);
//! Give one face of a body the support that `faceSpline` data describes, and grow or cut its
//! neighbours to meet it; see BRepSurfaceExt_FaceSupport.cxx.
Standard_EXPORT EvolutionData replaceFaceSupport(OcctKernel&                kernel,
                                                 uint32_t                   body,
                                                 uint32_t                   face,
                                                 const std::vector<double>& spline,
                                                 double                     tolerance,
                                                 const std::vector<int>&    hashes,
                                                 int                        bound);
//! The flat development of a plane, or of a rectangular cylindrical or conical face, in XY.
//! Fails on another surface or a cone that crosses its apex.
Standard_EXPORT EvolutionData unwrapFace(OcctKernel& kernel, uint32_t face, double tolerance);
//! The part of `wire` between the arc-length fractions `start` < `end` in [0, 1]. Fails where
//! the path has no length.
Standard_EXPORT EvolutionData trimPath(OcctKernel& kernel, uint32_t wire, double start, double end);
//! A face in the closed four-sided boundary of profile and guide `edges`, smooth to `firstFace`
//! and `lastFace` to the orders 0 to 2 (G0 to G2). Fails where the edges do not close.
Standard_EXPORT EvolutionData skinGuides(OcctKernel&                  kernel,
                                         const std::vector<uint32_t>& edges,
                                         uint32_t                     firstFace,
                                         uint32_t                     lastFace,
                                         int                          firstOrder,
                                         int                          lastOrder,
                                         double                       tolerance);
//! The parts of the section of bodies `first` and `second` that do not close: why a Boolean
//! found no region. Empty where every wire closes.
Standard_EXPORT EvolutionData openSections(OcctKernel& kernel,
                                           uint32_t    first,
                                           uint32_t    second,
                                           double      tolerance);
//! Each face of `body` moved by `distance` along its normal on its own, untrimmed: the shell an
//! offset would have to join.
Standard_EXPORT EvolutionData offsetShellFaces(OcctKernel& kernel, uint32_t body, double distance);
//! A solid made of `sheet` thickened by `thickness` (the sign picks the side). Fails where the
//! thickness passes a principal radius of curvature.
Standard_EXPORT EvolutionData thicken(OcctKernel&             kernel,
                                      uint32_t                sheet,
                                      double                  thickness,
                                      double                  tolerance,
                                      const std::vector<int>& hashes,
                                      int                     bound);
//! Why thickening `sheet` by `thickness` fails: always fails, naming the edges that shrink to
//! nothing or turn over (first order on curved faces).
Standard_EXPORT EvolutionData thickenFailure(OcctKernel& kernel, uint32_t sheet, double thickness);
//! The flat strip along a straight `reference` curve and its copy laid along the planar `target`
//! by arc length, wide enough for `moving`. Fails where the reference is not straight.
Standard_EXPORT EvolutionData deformSupports(OcctKernel& kernel,
                                             uint32_t    reference,
                                             uint32_t    target,
                                             uint32_t    moving,
                                             double      tolerance);
//! A fillet of `radius` on `edge` between its fractions `start` and `end`, along its tangent
//! chain with `propagate`; `law` optionally holds (fraction, radius) pairs from 0 to 1, flat.
Standard_EXPORT EvolutionData filletRange(OcctKernel&                kernel,
                                          uint32_t                   shape,
                                          uint32_t                   edge,
                                          double                     radius,
                                          double                     start,
                                          double                     end,
                                          bool                       propagate,
                                          const std::vector<int>&    hashes,
                                          int                        bound,
                                          const std::vector<double>& law = {});
//! A blend of `edge` between the planes square to it at the fractions `start` and `end` of it, the
//! rest of `shape` as it was: the blend of the edge alone (blending::Request: distances `first` and
//! `second` off `reference` or a chamfer's `angle` off it, a conic, G2 or chamfer section `law` of
//! `rho`, and where `sizes` holds (place, size) pairs, flat, a size that changes along it), limited
//! as filletRange's.
Standard_EXPORT EvolutionData blendRange(OcctKernel&                  kernel,
                                         uint32_t                     shape,
                                         uint32_t                     edge,
                                         double                       first,
                                         double                       second,
                                         blending::Law                law,
                                         double                       rho,
                                         uint32_t                     reference,
                                         double                       angle,
                                         double                       start,
                                         double                       end,
                                         const std::vector<int>&      hashes,
                                         int                          bound,
                                         const std::vector<double>&   sizes = {},
                                         std::vector<blending::Join>* joins = nullptr);
//! The indices, in the shape's map of edges from 0, of its free boundary edges.
Standard_EXPORT std::vector<int> boundaryEdgeIndices(OcctKernel& kernel, uint32_t shape);
//! Boundary `edges` of a planar face, or full circles of a cylindrical band, offset by `distance`
//! within the face. Fails on another surface or where the band collapses.
Standard_EXPORT EvolutionData offsetSheetEdge(OcctKernel&                  kernel,
                                              uint32_t                     face,
                                              const std::vector<uint32_t>& edges,
                                              double                       distance,
                                              double                       tolerance);
//! The face on the whole of its support; with `keepEdges`, divided along its old trim edges.
Standard_EXPORT EvolutionData untrimFace(OcctKernel& kernel,
                                         uint32_t    face,
                                         bool        keepEdges,
                                         double      tolerance);
//! Edges of constant u (`alongU`) or v at the `parameters`, fractions of the face's UV box,
//! clipped to the face. Fails where none meets it.
Standard_EXPORT EvolutionData isoparam(OcctKernel&                kernel,
                                       uint32_t                   face,
                                       bool                       alongU,
                                       const std::vector<double>& parameters,
                                       double                     tolerance);
//! An untrimmed rectangular sheet extended by `distance` at `side` with `continuity` 1 or 2.
//! Fails on a trimmed sheet or one with holes.
Standard_EXPORT EvolutionData extendSheet(OcctKernel&   kernel,
                                          uint32_t      face,
                                          SheetBoundary side,
                                          double        distance,
                                          int           continuity,
                                          double        tolerance);

//! An open wire extended at one end by `distance`, or up to `target` where one is given (0:
//! none). Fails on a closed wire.
Standard_EXPORT EvolutionData extendCurve(OcctKernel&    kernel,
                                          uint32_t       wire,
                                          bool           atStart,
                                          CurveExtension mode,
                                          double         distance,
                                          uint32_t       target,
                                          double         tolerance);

struct Vec3
{
  double x, y, z;
};

//! The sheet that `wire` sweeps along `direction`, and against it too with `both`.
Standard_EXPORT EvolutionData curveCurtain(OcctKernel& kernel,
                                           uint32_t    wire,
                                           Vec3        direction,
                                           bool        both);

//! `wire` projected onto the faces of `body` along `direction`, and against it too with `both`.
//! Fails where the projection does not cover the curve.
Standard_EXPORT EvolutionData
  projectCurve(OcctKernel& kernel, uint32_t wire, uint32_t body, Vec3 direction, bool both);

//! Every edge of `shape` projected along `normal` onto the plane through `origin`.
Standard_EXPORT EvolutionData projectPlane(OcctKernel& kernel,
                                           uint32_t    shape,
                                           Vec3        origin,
                                           Vec3        normal);

//! The silhouette of `shape` seen along `normal`, with `x` across the view: in space with `in3d`,
//! else laid on the plane through `origin`. Fails on a shape without faces.
Standard_EXPORT EvolutionData
  outline(OcctKernel& kernel, uint32_t shape, Vec3 origin, Vec3 normal, Vec3 x, bool in3d);

//! The closed outline at half `width` round an open planar `wire`, with round ends. Fails on a
//! closed or non-planar wire.
Standard_EXPORT EvolutionData slot(OcctKernel& kernel, uint32_t wire, double width, Vec3 normal);

//! A wire winding `turns` times about the axis from `start`'s side: `radius` away at first and
//! `pitch` further each turn, rising `height` in all (BRepCurveExt_Spiral.cxx). Approximated to
//! `tolerance`.
Standard_EXPORT EvolutionData spiral(OcctKernel& kernel,
                                     Vec3        origin,
                                     Vec3        axis,
                                     Vec3        start,
                                     double      radius,
                                     double      height,
                                     double      turns,
                                     double      pitch,
                                     bool        leftHanded,
                                     double      tolerance);

//! The least distance between two shapes and the nearest point on each: d, then x, y, z twice.
Standard_EXPORT std::vector<double> minimumDistance(OcctKernel& kernel,
                                                    uint32_t    first,
                                                    uint32_t    second);
//! At `samples` places along `edge`, shared by two faces: x, y, z, gap, normal angle (degrees),
//! curvature difference, relative difference. Fails where the edge is not on both faces.
Standard_EXPORT std::vector<double> faceContinuity(OcctKernel& kernel,
                                                   uint32_t    first,
                                                   uint32_t    second,
                                                   uint32_t    edge,
                                                   int         samples,
                                                   uint32_t    otherEdge = 0);

//! An unmodified primitive's PrimitiveKind, origin, x, y and z axes, and dimensions a, b, c
//! (box sides; or diameter and height; or diameter). Empty for another shape.
Standard_EXPORT std::vector<double> primitiveDimensions(OcctKernel& kernel,
                                                        uint32_t    shape,
                                                        double      tolerance);
//! The primitive of `kind` rebuilt with dimensions `a`, `b`, `c` in its own frame. Fails where
//! `shape` is not that unmodified primitive or a dimension is not positive.
Standard_EXPORT EvolutionData setPrimitiveDimensions(OcctKernel&             kernel,
                                                     uint32_t                shape,
                                                     PrimitiveKind           kind,
                                                     double                  a,
                                                     double                  b,
                                                     double                  c,
                                                     double                  tolerance,
                                                     const std::vector<int>& hashes,
                                                     int                     bound);

//! A circle in the plane tangent to `first` and `second` with `radius`, or to three edges, the
//! candidate nearest `near`. Fails where none is found.
Standard_EXPORT EvolutionData tangentCircle(OcctKernel& kernel,
                                            uint32_t    first,
                                            uint32_t    second,
                                            uint32_t    third,
                                            double      radius,
                                            Vec3        origin,
                                            Vec3        normal,
                                            Vec3        x,
                                            Vec3        near,
                                            double      tolerance);

//! `wire` refitted as splines of `degree` 1 to 9, within `tolerance` (by edge with `keepCorners`)
//! or with `count` poles. Fails where the refit misses the tolerance.
Standard_EXPORT EvolutionData rebuildCurve(OcctKernel&  kernel,
                                           uint32_t     wire,
                                           CurveRebuild mode,
                                           double       tolerance,
                                           bool         keepCorners,
                                           int          degree,
                                           int          count);

//! A chamfer of `edges`: by `distance1`, or off `reference` by two distances or a distance and
//! `angle` (radians). Fails on distances that are not positive or tangent edges.
Standard_EXPORT EvolutionData chamfer(OcctKernel&                  kernel,
                                      uint32_t                     shape,
                                      const std::vector<uint32_t>& edges,
                                      double                       distance1,
                                      double                       distance2,
                                      uint32_t                     reference,
                                      const std::vector<int>&      hashes,
                                      int                          bound,
                                      double                       angle,
                                      bool                         propagate);

//! For each blend face of `shape`: its face index, kind, neighbour count and the indices of its
//! blend neighbours across edges.
Standard_EXPORT std::vector<double> blendTopology(OcctKernel& kernel,
                                                  uint32_t    shape,
                                                  double      tolerance);
//! A cylindrical or toroidal fillet face's radius, centre and the face's middle point (x, y, z
//! each). Empty for another face.
Standard_EXPORT std::vector<double> filletDimension(OcctKernel& kernel,
                                                    uint32_t    shape,
                                                    uint32_t    face,
                                                    double      tolerance);
//! The fillet `face` replaced by one of `radius`. Fails where the face is not a fillet or its
//! removal does not leave one edge.
Standard_EXPORT EvolutionData refillet(OcctKernel&             kernel,
                                       uint32_t                shape,
                                       uint32_t                face,
                                       double                  radius,
                                       double                  tolerance,
                                       const std::vector<int>& hashes,
                                       int                     bound);

//! At each arc-length fraction of `wire`: point, tangent, normal and binormal; the frame
//! rotation-minimising with `transport`. Fails at a cusp.
Standard_EXPORT std::vector<double> curveLocations(OcctKernel&                kernel,
                                                   uint32_t                   wire,
                                                   const std::vector<double>& fractions,
                                                   bool                       transport);

//! What of `shape` lies on the side of the plane that `normal` points to (`positive`), or on the
//! other. Fails on a zero normal.
Standard_EXPORT EvolutionData clipPlane(OcctKernel&             kernel,
                                        uint32_t                shape,
                                        Vec3                    origin,
                                        Vec3                    normal,
                                        bool                    positive,
                                        const std::vector<int>& hashes,
                                        int                     bound);

//! A planar region `face` grown (positive `distance`) or shrunk. Fails where nothing of it is
//! left.
Standard_EXPORT EvolutionData
  offsetRegion(OcctKernel& kernel, uint32_t face, double distance, OffsetGap gap, double tolerance);
//! A planar `wire` offset by `distance`, outwards or to the right where positive.
Standard_EXPORT EvolutionData offsetCurve(OcctKernel& kernel,
                                          uint32_t    wire,
                                          double      distance,
                                          OffsetGap   gap);
//! offsetCurve, with the right of an open wire taken against `up`.
Standard_EXPORT EvolutionData
  offsetCurveOn(OcctKernel& kernel, uint32_t wire, double distance, OffsetGap gap, Vec3 up);

//! `wire` with new vertices `distance` along it on either side of `vertex`. Fails where the
//! distance runs past an end.
Standard_EXPORT EvolutionData offsetVertex(OcctKernel& kernel,
                                           uint32_t    wire,
                                           uint32_t    vertex,
                                           double      distance,
                                           double      tolerance);

//! An arc in the plane of `normal` leaving `start` along `tangent` and ending tangent to `target`,
//! the one nearest `near`. Fails where none exists.
Standard_EXPORT EvolutionData tangentArcToCurve(OcctKernel& kernel,
                                                Vec3        start,
                                                Vec3        tangent,
                                                uint32_t    target,
                                                Vec3        normal,
                                                Vec3        near,
                                                double      tolerance);

//! Push-pull: move `faces` of a solid or a sheet along their outward normals by
//! `distance` (negative moves into the material). Neighbouring faces are
//! extended or trimmed to meet them; a blend tangent to a flat face moves with it, and a
//! move may remove faces or open the body.
Standard_EXPORT EvolutionData offsetFaces(OcctKernel&                  kernel,
                                          uint32_t                     body,
                                          const std::vector<uint32_t>& faces,
                                          double                       distance,
                                          double                       tolerance,
                                          const std::vector<int>&      hashes,
                                          int                          bound,
                                          FaceGrow                     grow = FaceGrow::None);

//! Taper `faces` by `angle` (radians) about their intersection with the
//! neutral plane
//! (`planeOrigin`, `planeNormal`). Material is removed on the `pullDirection`
//! side of the plane for a positive angle and added for a negative one.
Standard_EXPORT EvolutionData draftFaces(OcctKernel&                  kernel,
                                         uint32_t                     shape,
                                         const std::vector<uint32_t>& faces,
                                         double                       angle,
                                         Vec3                         pullDirection,
                                         Vec3                         planeOrigin,
                                         Vec3                         planeNormal,
                                         const std::vector<int>&      hashes,
                                         int                          bound);

//! Remove `faces` from a solid and close the gap by extending their
//! neighbours.
Standard_EXPORT EvolutionData defeature(OcctKernel&                  kernel,
                                        uint32_t                     shape,
                                        const std::vector<uint32_t>& faces,
                                        const std::vector<int>&      hashes,
                                        int                          bound);

//! Fillet with a radius that varies along each edge. Edge i takes `lawSizes[i]`
//! consecutive (parameter, radius) pairs from `laws`; the parameter runs from 0
//! at the start of the edge's curve to 1 at its end.
Standard_EXPORT EvolutionData filletVariable(OcctKernel&                  kernel,
                                             uint32_t                     shape,
                                             const std::vector<uint32_t>& edges,
                                             const std::vector<int>&      lawSizes,
                                             const std::vector<double>&   laws,
                                             const std::vector<int>&      hashes,
                                             int                          bound);

//! Scale each world axis about a fixed center, reporting the actual OCCT face
//! history.
Standard_EXPORT EvolutionData scaleAxes(OcctKernel&             kernel,
                                        uint32_t                shape,
                                        Vec3                    factors,
                                        Vec3                    center,
                                        const std::vector<int>& hashes,
                                        int                     bound);

//! All bounded cells of coplanar curves, including intersections and nested
//! loops.
Standard_EXPORT EvolutionData regions(OcctKernel&                  kernel,
                                      const std::vector<uint32_t>& curves,
                                      Vec3                         origin,
                                      Vec3                         normal,
                                      double                       tolerance);

//! Signed volume using adaptive integration of curved faces.
Standard_EXPORT double volume(OcctKernel& kernel, uint32_t shape);
//! Volume by OpenCascade's fixed quadrature, without an error bound.
Standard_EXPORT double volumeEstimate(OcctKernel& kernel, uint32_t shape);
//! Area by OpenCascade's fixed quadrature, without an error bound.
Standard_EXPORT double surfaceAreaEstimate(OcctKernel& kernel, uint32_t shape);
//! A box of the exact geometry: xmin, ymin, zmin, xmax, ymax, zmax. Fails on a shape without
//! geometry.
Standard_EXPORT std::vector<double> boundingBox(OcctKernel& kernel, uint32_t shape);
//! What the face separation pre-pass proves of a shape (BRepSurfaceExt_FaceSeparation.hxx): its
//! face pairs (each face with itself included), those proven, its faces, those a check still needs,
//! and the parts it checks them in.
Standard_EXPORT std::vector<double> faceSeparation(OcctKernel& kernel, uint32_t shape);
//! The face separation pre-pass checked: its proven pairs, those that do meet, aborted checks,
//! and whether the shape crosses itself with the pre-pass and without (1 or 0 each).
Standard_EXPORT std::vector<double> faceSeparationAudit(OcctKernel& kernel, uint32_t shape);
//! The signed volume by adaptive integration, then 1 where every face converged, else 0.
Standard_EXPORT std::vector<double> volumeResult(OcctKernel& kernel, uint32_t shape);
//! The area by adaptive integration, then 1 where every face converged, else 0.
Standard_EXPORT std::vector<double> surfaceAreaResult(OcctKernel& kernel, uint32_t shape);
//! The area centroid: the facade's, except of faces of a support other than a plane, cylinder,
//! cone, sphere or torus trimmed to less than their parameter rectangle, which are integrated
//! over their trim only.
Standard_EXPORT std::vector<double> surfaceCenter(OcctKernel& kernel, uint32_t shape);
//! The area by adaptive integration.
Standard_EXPORT double surfaceArea(OcctKernel& kernel, uint32_t shape);

//! Position, unit tangent and curvature vector, in groups of nine, from
//! topological edge start.
Standard_EXPORT std::vector<double> curvature(OcctKernel&                kernel,
                                              uint32_t                   edge,
                                              const std::vector<double>& parameters);

//! Hollow the entire solid without opening its outer skin.
Standard_EXPORT EvolutionData hollow(OcctKernel&             kernel,
                                     uint32_t                shape,
                                     double                  thickness,
                                     double                  tolerance,
                                     const std::vector<int>& hashes,
                                     int                     bound);

//! A shape-preserving edit of each edge of `wire`: its degree raised by `amount` (at most 24),
//! its poles doubled, or redundant edges merged within `tolerance`.
Standard_EXPORT EvolutionData refineCurve(OcctKernel&     kernel,
                                          uint32_t        wire,
                                          CurveRefinement mode,
                                          int             amount,
                                          double          tolerance);

//! A B-spline face through a `rows` by `columns` grid of `coordinates` (x, y, z in turn, by row).
//! Fails where it misses `tolerance`.
Standard_EXPORT EvolutionData surfaceThroughGrid(OcctKernel&                kernel,
                                                 const std::vector<double>& coordinates,
                                                 int                        rows,
                                                 int                        columns,
                                                 double                     tolerance);
//! A face filling the closed `wire`. Fails where it misses the boundary by more than `tolerance`.
Standard_EXPORT EvolutionData fillBoundary(OcctKernel& kernel, uint32_t wire, double tolerance);

//! What `patch` built and how it measures: the construction (0 a plane, 1 the Coons patch of
//! four sides, 2 four sides fitted, 3 a trimmed fitted surface, 4 OpenCascade's plate), the
//! degrees and the pole counts in u and v; then per boundary edge the largest gap, angle in
//! degrees, share and difference of curvature; then the largest gap per guide and per point.
struct PatchResult
{
  EvolutionData       evolution;
  std::vector<double> values;
};

//! How `patch` constructs the face.
enum class PatchMethod
{
  Auto   = 0, //!< its own choice
  Fitted = 1, //!< never the plane's shortcut
  Plate  = 2  //!< OpenCascade's plate
};

//! Whether the face is made of four sides untrimmed.
enum class PatchQuad
{
  Auto = -1,
  Off  = 0,
  On   = 1
};

//! The flow of the face's parameters at its sides.
enum class PatchFlow
{
  Natural  = 0,
  Normal   = 1,
  Next     = 2,
  Adjacent = 3
};

//! The least number of fitted spans: Auto none, High 10, Max 16.
enum class PatchQuality
{
  Auto = 0,
  High = 1,
  Max  = 2
};

//! Controls on this fill; the defaults keep Patch's established automatic construction.
struct PatchControls
{
  PatchQuad    quad    = PatchQuad::Auto;
  PatchFlow    flow    = PatchFlow::Natural;
  PatchQuality quality = PatchQuality::Auto;
  bool         trim    = true;
  bool         open = false; // Edges that do not close are chained, the gaps closed straight in UV.
};

//! One face inside the closed loop of `edges`. `supports[i]` is the face that edge i is to be
//! smooth with to `orders[i]` (0, 1 or 2 for G0 to G2), or 0 for none; `guides` are edges and
//! `points` coordinates in threes that the face is to pass through. Fails where what it built
//! is measured to miss `tolerance`, `g1Degrees` or `g2Share`.
Standard_EXPORT PatchResult patch(OcctKernel&                  kernel,
                                  const std::vector<uint32_t>& edges,
                                  const std::vector<uint32_t>& supports,
                                  const std::vector<int>&      orders,
                                  const std::vector<uint32_t>& guides,
                                  const std::vector<double>&   points,
                                  double                       tolerance,
                                  double                       g1Degrees,
                                  double                       g2Share,
                                  PatchMethod                  method,
                                  double                       flatness = -1,
                                  const PatchControls&         controls = {});

//! How fair a face is, sampled `samples` times each way over its parameters: its area; the
//! integral of the squared principal curvatures; the integral of the squared change of the
//! principal curvatures per length, times the area; and the largest curvature. The second and
//! third are plain numbers, the same at every scale.
Standard_EXPORT std::vector<double> faceFairness(OcctKernel& kernel, uint32_t face, int samples);

//! Deform Curve in UVN: see BRepCurveExt_EditMore.cxx.
Standard_EXPORT EvolutionData deformCurveUVN(OcctKernel& kernel,
                                             uint32_t    wire,
                                             uint32_t    reference,
                                             uint32_t    target,
                                             Vec3        scale,
                                             Vec3        offset,
                                             bool        swapUV,
                                             bool        mirror,
                                             bool        flipNormal,
                                             double      tolerance);

//! Rebuild an open curve as one spline that is near it and smooth, by `weight` (above 0, at
//! most 1: 1 is nearest): see BRepCurveExt_EditMore.cxx.
Standard_EXPORT EvolutionData
  rebuildCurveWeighted(OcctKernel& kernel, uint32_t wire, int degree, int count, double weight);

//! Give each spline edge that uses a part of its curve a spline of exactly that part.
Standard_EXPORT EvolutionData untrimCurve(OcctKernel& kernel, uint32_t wire);

//! `projectCurve` with what is hidden from the curve left out.
Standard_EXPORT EvolutionData
  projectCurveVisible(OcctKernel& kernel, uint32_t wire, uint32_t body, Vec3 direction, bool both);

//! A sweep along a smooth open path with rotation-minimising frames, a twist (radians) and a
//! scale at the far end; or, with one to three guides and neither twist nor scale, the section
//! placed by the guides; `guideMethod` matters for two guides. See BRepSkinExt_SweepFramed.cxx.
Standard_EXPORT EvolutionData sweepFramed(OcctKernel&                  kernel,
                                          uint32_t                     profile,
                                          uint32_t                     path,
                                          double                       twist,
                                          double                       scale,
                                          const std::vector<uint32_t>& guides,
                                          GuideMethod                  guideMethod,
                                          double                       tolerance,
                                          bool                         solid);

//! A skin through `sections`, from and to a point vertex where given (0: none): smooth through
//! three or more, else ruled. Fails where the skin is not valid. See BRepSkinExt_Gordon.cxx.
Standard_EXPORT EvolutionData skinThrough(OcctKernel&                  kernel,
                                          const std::vector<uint32_t>& sections,
                                          uint32_t                     startVertex,
                                          uint32_t                     endVertex,
                                          bool                         solid,
                                          bool                         ruled);
//! A skin through two or more profiles along two or more guides, each an edge or a wire, as a
//! Gordon surface; with `solid`, closed flat end profiles are capped. See BRepSkinExt_Gordon.cxx.
Standard_EXPORT EvolutionData skinGordon(OcctKernel&                  kernel,
                                         const std::vector<uint32_t>& profiles,
                                         const std::vector<uint32_t>& guides,
                                         bool                         solid,
                                         double                       tolerance);

//! A skin through open profiles along guides that cross them anywhere, each guide's move off the
//! profiles' skin blended to the next by a cubic without slope at either. See
//! BRepSkinExt_Gordon.cxx.
Standard_EXPORT EvolutionData skinGuided(OcctKernel&                  kernel,
                                         const std::vector<uint32_t>& profiles,
                                         const std::vector<uint32_t>& guides,
                                         double                       tolerance);

//! Edges on a face of a body along lines of constant u or v: see BRepBooleanExt_ImprintMore.cxx.
Standard_EXPORT EvolutionData isoparamEdges(OcctKernel&                kernel,
                                            uint32_t                   shape,
                                            uint32_t                   face,
                                            bool                       alongU,
                                            const std::vector<double>& parameters,
                                            bool                       subdivide,
                                            double                     tolerance,
                                            const std::vector<int>&    hashes,
                                            int                        bound);

//! Imprint curves on a body by a direction or along the faces' normals, optionally without
//! what is hidden, and with loose ends carried on to the faces' edges: see
//! BRepBooleanExt_ImprintMore.cxx.
Standard_EXPORT EvolutionData imprintCurves(OcctKernel&                  kernel,
                                            uint32_t                     shape,
                                            const std::vector<uint32_t>& curves,
                                            ImprintMethod                method,
                                            Vec3                         direction,
                                            bool                         both,
                                            bool                         hideOccluded,
                                            ImprintCompletion            completion,
                                            const std::vector<int>&      hashes,
                                            int                          bound);

//! Complete Edge: the loose ends of edges inside faces run on across their faces until they
//! meet an edge, and the body is split along them: see BRepBooleanExt_ImprintMore.cxx.
Standard_EXPORT EvolutionData completeEdges(OcctKernel&                  kernel,
                                            uint32_t                     shape,
                                            const std::vector<uint32_t>& edges,
                                            double                       tolerance,
                                            const std::vector<int>&      hashes,
                                            int                          bound);

//! Every face of a solid moved by `distance` along its outward normal (negative: inwards).
//! Faces that the distance shrinks to nothing are removed first; a face that is too tight in
//! part only is refused, naming it.
Standard_EXPORT EvolutionData offsetSolid(OcctKernel&             kernel,
                                          uint32_t                shape,
                                          double                  distance,
                                          double                  tolerance,
                                          OffsetJoin              join,
                                          const std::vector<int>& hashes,
                                          int                     bound);

//! The faces of a body split along curves that lie on them, and only those faces (the Knife):
//! see BRepBooleanExt_ImprintMore.cxx.
Standard_EXPORT EvolutionData splitOnFaces(OcctKernel&                  kernel,
                                           uint32_t                     shape,
                                           const std::vector<uint32_t>& curves,
                                           double                       tolerance,
                                           const std::vector<int>&      hashes,
                                           int                          bound);

//! Spatial radius stations and directed limits along selected tangent contours.
Standard_EXPORT EvolutionData filletStations(OcctKernel&                  kernel,
                                             uint32_t                     shape,
                                             const std::vector<uint32_t>& edges,
                                             double                       radius,
                                             const std::vector<uint32_t>& stationEdges,
                                             const std::vector<double>&   stations,
                                             const std::vector<uint32_t>& limitEdges,
                                             const std::vector<double>&   limits,
                                             const std::vector<int>&      hashes,
                                             int                          bound);

//! A planar `profile` swept along `path`, scaled evenly to `scale` at its end; with `parallel`,
//! the sections stay parallel to the profile. Fails on a profile that is not planar.
Standard_EXPORT EvolutionData sweepSection(OcctKernel& kernel,
                                           uint32_t    profile,
                                           uint32_t    path,
                                           double      scale,
                                           SweepCorner corner,
                                           bool        parallel,
                                           bool        solid,
                                           double      tolerance);
//! One `face` of a body divided by `curves` clipped to it. Fails where they lie outside the face
//! or do not divide it.
Standard_EXPORT EvolutionData splitFaceEdges(OcctKernel&                  kernel,
                                             uint32_t                     shape,
                                             uint32_t                     face,
                                             const std::vector<uint32_t>& curves,
                                             double                       tolerance,
                                             const std::vector<int>&      hashes,
                                             int                          bound);

} // namespace BRepExt
