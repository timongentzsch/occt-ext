// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepToolsExt_Session.hxx>

// The operations this package adds to the kernel. Each takes and returns shape
// handles, reports face history like upstream's *WithHistory methods (see
// BRepToolsExt_History.hxx) and fails through BRepExt::fail (see GeomAlgoExt_Failure.hxx).

namespace BRepExt
{
namespace blending { struct Join; }
struct Failure;
Standard_EXPORT void          clearTriangulation(OcctKernel&, uint32_t);
Standard_EXPORT EvolutionData makeBox(OcctKernel&, double, double, double);
Standard_EXPORT EvolutionData booleanChecked(OcctKernel&,
                                             int,
                                             const std::vector<uint32_t>&,
                                             const std::vector<uint32_t>&,
                                             int,
                                             double,
                                             double,
                                             const std::vector<int>&,
                                             int);
//! A selection of the cells of one General Fuse (BRepBooleanExt_Cells.cxx: Select): 0 union, 1
//! difference, 2 intersection, 3 slice, 4 imprint, 5 what is inside every operand.
Standard_EXPORT EvolutionData booleanCells(OcctKernel&,
                                           int                          select,
                                           const std::vector<uint32_t>& sources,
                                           const std::vector<uint32_t>& tools,
                                           double                       fuzzy,
                                           const std::vector<int>&      hashes,
                                           int                          bound);
struct BooleanAttempt;
//! What each rung of the ladder did for that selection, whether or not one was accepted.
Standard_EXPORT std::vector<BooleanAttempt> booleanReport(OcctKernel&,
                                                          int                          select,
                                                          const std::vector<uint32_t>& sources,
                                                          const std::vector<uint32_t>& tools,
                                                          double                       fuzzy);
Standard_EXPORT EvolutionData     shell(OcctKernel&,
                                        uint32_t,
                                        const std::vector<uint32_t>&,
                                        double,
                                        double,
                                        int,
                                        const std::vector<int>&,
                                        int,
                                        bool);

Standard_EXPORT std::vector<int> surfaceSplineInfo(OcctKernel&, uint32_t);
struct Vec3;
Standard_EXPORT EvolutionData interpolateClosed(OcctKernel&, const std::vector<double>&);
Standard_EXPORT EvolutionData skinPeriodic(OcctKernel&, const std::vector<uint32_t>&, double);
Standard_EXPORT EvolutionData importIges(OcctKernel&, const std::string&);
Standard_EXPORT std::string   exportIges(OcctKernel&, uint32_t);
Standard_EXPORT EvolutionData rebuildFace(OcctKernel&, uint32_t, int, int, int, double);
//! Rebuild one face of a body in place: its support is refitted on the face's UV box and its
//! edges stay shared with the neighbours, within the fit tolerance.
Standard_EXPORT EvolutionData rebuildBodyFace(OcctKernel&,
                                              uint32_t body,
                                              uint32_t face,
                                              int,
                                              int,
                                              int,
                                              double,
                                              const std::vector<int>& hashes,
                                              int                     bound);
//! A B-spline or Bezier face support as exact world-space B-spline data: degrees,
//! knot and pole counts, knots with
//! multiplicities in U then V, then each pole as x, y, z, weight, U-major.
Standard_EXPORT std::vector<double> faceSpline(OcctKernel&, uint32_t face);
//! A face on the whole of the spline that `faceSpline` describes.
Standard_EXPORT EvolutionData splineFace(OcctKernel&, const std::vector<double>&);
//! Give one face of a body the support that `faceSpline` data describes, and grow or cut its
//! neighbours to meet it; see BRepSurfaceExt_FaceSupport.cxx.
Standard_EXPORT EvolutionData replaceFaceSupport(OcctKernel&,
                                                 uint32_t                   body,
                                                 uint32_t                   face,
                                                 const std::vector<double>& spline,
                                                 double                     tolerance,
                                                 const std::vector<int>&    hashes,
                                                 int                        bound);
Standard_EXPORT EvolutionData unwrapFace(OcctKernel&, uint32_t, double);
Standard_EXPORT EvolutionData trimPath(OcctKernel&, uint32_t, double, double);
Standard_EXPORT EvolutionData skinGuides(OcctKernel&,
                                         const std::vector<uint32_t>& edges,
                                         uint32_t                     firstFace,
                                         uint32_t                     lastFace,
                                         int                          firstOrder,
                                         int                          lastOrder,
                                         double                       tolerance);
Standard_EXPORT EvolutionData openSections(OcctKernel&, uint32_t a, uint32_t b, double tolerance);
Standard_EXPORT EvolutionData offsetShellFaces(OcctKernel&, uint32_t body, double distance);
Standard_EXPORT EvolutionData
  thicken(OcctKernel&, uint32_t, double, double, const std::vector<int>&, int);
Standard_EXPORT EvolutionData thickenFailure(OcctKernel&, uint32_t sheet, double thickness);
Standard_EXPORT EvolutionData deformSupports(OcctKernel&,
                                             uint32_t reference,
                                             uint32_t target,
                                             uint32_t moving,
                                             double   tolerance);
Standard_EXPORT EvolutionData filletRange(OcctKernel&,
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
Standard_EXPORT EvolutionData blendRange(OcctKernel&,
                                         uint32_t                   shape,
                                         uint32_t                   edge,
                                         double                     first,
                                         double                     second,
                                         int                        law,
                                         double                     rho,
                                         uint32_t                   reference,
                                         double                     angle,
                                         double                     start,
                                         double                     end,
                                         const std::vector<int>&    hashes,
                                         int                        bound,
                                         const std::vector<double>& sizes = {},
                                         std::vector<blending::Join>* joins = nullptr);
Standard_EXPORT std::vector<int> boundaryEdgeIndices(OcctKernel&, uint32_t shape);
Standard_EXPORT EvolutionData    offsetSheetEdge(OcctKernel&,
                                                 uint32_t                     face,
                                                 const std::vector<uint32_t>& edges,
                                                 double                       distance,
                                                 double                       tolerance);
Standard_EXPORT EvolutionData    untrimFace(OcctKernel&,
                                            uint32_t face,
                                            bool     keepEdges,
                                            double   tolerance);
Standard_EXPORT EvolutionData    isoparam(OcctKernel&,
                                          uint32_t                   face,
                                          bool                       alongU,
                                          const std::vector<double>& parameters,
                                          double                     tolerance);
Standard_EXPORT EvolutionData    extendSheet(OcctKernel&,
                                             uint32_t face,
                                             int      side,
                                             double   distance,
                                             int      continuity,
                                             double   tolerance);

Standard_EXPORT EvolutionData extendCurve(OcctKernel& kernel,
                                          uint32_t    wire,
                                          bool        atStart,
                                          int         mode,
                                          double      distance,
                                          uint32_t    target,
                                          double      tolerance);

struct Vec3
{
  double x, y, z;
};

Standard_EXPORT EvolutionData curveCurtain(OcctKernel& kernel,
                                           uint32_t    wire,
                                           Vec3        direction,
                                           bool        both);

Standard_EXPORT EvolutionData
  projectCurve(OcctKernel& kernel, uint32_t wire, uint32_t body, Vec3 direction, bool both);

Standard_EXPORT EvolutionData projectPlane(OcctKernel& kernel,
                                           uint32_t    shape,
                                           Vec3        origin,
                                           Vec3        normal);

Standard_EXPORT EvolutionData
  outline(OcctKernel& kernel, uint32_t shape, Vec3 origin, Vec3 normal, Vec3 x, bool in3d);

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

Standard_EXPORT std::vector<double> minimumDistance(OcctKernel& kernel,
                                                    uint32_t    first,
                                                    uint32_t    second);
Standard_EXPORT std::vector<double> faceContinuity(OcctKernel& kernel,
                                                   uint32_t    first,
                                                   uint32_t    second,
                                                   uint32_t    edge,
                                                   int         samples,
                                                   uint32_t    otherEdge = 0);

Standard_EXPORT std::vector<double> primitiveDimensions(OcctKernel& kernel,
                                                        uint32_t    shape,
                                                        double      tolerance);
Standard_EXPORT EvolutionData       setPrimitiveDimensions(OcctKernel&             kernel,
                                                           uint32_t                shape,
                                                           int                     kind,
                                                           double                  a,
                                                           double                  b,
                                                           double                  c,
                                                           double                  tolerance,
                                                           const std::vector<int>& hashes,
                                                           int                     bound);

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

Standard_EXPORT EvolutionData rebuildCurve(OcctKernel& kernel,
                                           uint32_t    wire,
                                           int         mode,
                                           double      tolerance,
                                           bool        keepCorners,
                                           int         degree,
                                           int         count);

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

Standard_EXPORT std::vector<double> blendTopology(OcctKernel& kernel,
                                                  uint32_t    shape,
                                                  double      tolerance);
Standard_EXPORT std::vector<double> filletDimension(OcctKernel& kernel,
                                                    uint32_t    shape,
                                                    uint32_t    face,
                                                    double      tolerance);
Standard_EXPORT EvolutionData       refillet(OcctKernel&             kernel,
                                             uint32_t                shape,
                                             uint32_t                face,
                                             double                  radius,
                                             double                  tolerance,
                                             const std::vector<int>& hashes,
                                             int                     bound);

Standard_EXPORT std::vector<double> curveLocations(OcctKernel&                kernel,
                                                   uint32_t                   wire,
                                                   const std::vector<double>& fractions,
                                                   bool                       transport);

Standard_EXPORT EvolutionData clipPlane(OcctKernel&             kernel,
                                        uint32_t                shape,
                                        Vec3                    origin,
                                        Vec3                    normal,
                                        bool                    positive,
                                        const std::vector<int>& hashes,
                                        int                     bound);

Standard_EXPORT EvolutionData
  offsetRegion(OcctKernel& kernel, uint32_t face, double distance, int gap, double tolerance);
Standard_EXPORT EvolutionData offsetCurve(OcctKernel& kernel,
                                          uint32_t    wire,
                                          double      distance,
                                          int         gap);
Standard_EXPORT EvolutionData
  offsetCurveOn(OcctKernel& kernel, uint32_t wire, double distance, int gap, Vec3 up);

Standard_EXPORT EvolutionData offsetVertex(OcctKernel& kernel,
                                           uint32_t    wire,
                                           uint32_t    vertex,
                                           double      distance,
                                           double      tolerance);

Standard_EXPORT EvolutionData tangentArcToCurve(OcctKernel& kernel,
                                                Vec3        start,
                                                Vec3        tangent,
                                                uint32_t    target,
                                                Vec3        normal,
                                                Vec3        near,
                                                double      tolerance);

//! Push-pull: move `faceIds` of a solid or a sheet along their outward normals by
//! `distance` (negative moves into the material). Neighbouring faces are
//! extended or trimmed to meet them; a blend tangent to a flat face moves with it, and a
//! move may remove faces or open the body. `grow`, for one flat face pulled out to a wall
//! that faces the same way: 0 goes on past it, 1 stops flush with it, 2 goes on with it.
Standard_EXPORT EvolutionData offsetFaces(OcctKernel&                  kernel,
                                          uint32_t                     solidId,
                                          const std::vector<uint32_t>& faceIds,
                                          double                       distance,
                                          double                       tolerance,
                                          const std::vector<int>&      inputFaceHashes,
                                          int                          hashUpperBound,
                                          int                          grow = 0);

//! Taper `faceIds` by `angle` (radians) about their intersection with the
//! neutral plane
//! (`planeOrigin`, `planeNormal`). Material is removed on the `pullDirection`
//! side of the plane for a positive angle and added for a negative one.
Standard_EXPORT EvolutionData draftFaces(OcctKernel&                  kernel,
                                         uint32_t                     shapeId,
                                         const std::vector<uint32_t>& faceIds,
                                         double                       angle,
                                         Vec3                         pullDirection,
                                         Vec3                         planeOrigin,
                                         Vec3                         planeNormal,
                                         const std::vector<int>&      inputFaceHashes,
                                         int                          hashUpperBound);

//! Remove `faceIds` from a solid and close the gap by extending their
//! neighbours.
Standard_EXPORT EvolutionData defeature(OcctKernel&                  kernel,
                                        uint32_t                     shapeId,
                                        const std::vector<uint32_t>& faceIds,
                                        const std::vector<int>&      inputFaceHashes,
                                        int                          hashUpperBound);

//! Fillet with a radius that varies along each edge. Edge i takes `lawSizes[i]`
//! consecutive (parameter, radius) pairs from `laws`; the parameter runs from 0
//! at the start of the edge's curve to 1 at its end.
Standard_EXPORT EvolutionData filletVariable(OcctKernel&                  kernel,
                                             uint32_t                     shapeId,
                                             const std::vector<uint32_t>& edgeIds,
                                             const std::vector<int>&      lawSizes,
                                             const std::vector<double>&   laws,
                                             const std::vector<int>&      inputFaceHashes,
                                             int                          hashUpperBound);

//! Scale each world axis about a fixed center, reporting the actual OCCT face
//! history.
Standard_EXPORT EvolutionData scaleAxes(OcctKernel&             kernel,
                                        uint32_t                shapeId,
                                        Vec3                    factors,
                                        Vec3                    center,
                                        const std::vector<int>& inputFaceHashes,
                                        int                     hashUpperBound);

//! All bounded cells of coplanar curves, including intersections and nested
//! loops.
Standard_EXPORT EvolutionData regions(OcctKernel&                  kernel,
                                      const std::vector<uint32_t>& curveIds,
                                      Vec3                         origin,
                                      Vec3                         normal,
                                      double                       tolerance);

//! Signed volume using adaptive integration of curved faces.
Standard_EXPORT double volume(OcctKernel& kernel, uint32_t shapeId);
Standard_EXPORT double volumeEstimate(OcctKernel& kernel, uint32_t shapeId);
Standard_EXPORT double surfaceAreaEstimate(OcctKernel& kernel, uint32_t shapeId);
Standard_EXPORT std::vector<double> boundingBox(OcctKernel&, uint32_t);
//! What the face separation pre-pass proves of a shape (BRepSurfaceExt_FaceSeparation.hxx): its
//! face pairs (each face with itself included), those proven, its faces, those a check still needs,
//! and the parts it checks them in.
Standard_EXPORT std::vector<double> faceSeparation(OcctKernel&, uint32_t);
Standard_EXPORT std::vector<double> faceSeparationAudit(OcctKernel&, uint32_t);
Standard_EXPORT std::vector<double> volumeResult(OcctKernel&, uint32_t);
Standard_EXPORT std::vector<double> surfaceAreaResult(OcctKernel&, uint32_t);
//! The area centroid: the facade's, except of faces of a support other than a plane, cylinder,
//! cone, sphere or torus trimmed to less than their parameter rectangle, which are integrated
//! over their trim only.
Standard_EXPORT std::vector<double> surfaceCenter(OcctKernel&, uint32_t);
Standard_EXPORT double              surfaceArea(OcctKernel& kernel, uint32_t shapeId);

//! Position, unit tangent and curvature vector, in groups of nine, from
//! topological edge start.
Standard_EXPORT std::vector<double> curvature(OcctKernel&                kernel,
                                              uint32_t                   edgeId,
                                              const std::vector<double>& parameters);

//! Hollow the entire solid without opening its outer skin.
Standard_EXPORT EvolutionData hollow(OcctKernel&             kernel,
                                     uint32_t                shapeId,
                                     double                  thickness,
                                     double                  tolerance,
                                     const std::vector<int>& inputFaceHashes,
                                     int                     hashUpperBound);

//! Shape-preserving curve edits: 0 raises degree, 1 doubles CVs, 2 removes
//! redundant topology.
Standard_EXPORT EvolutionData
  refineCurve(OcctKernel& kernel, uint32_t wireId, int mode, int amount, double tolerance);

Standard_EXPORT EvolutionData surfaceThroughGrid(OcctKernel&                kernel,
                                                 const std::vector<double>& coordinates,
                                                 int                        rows,
                                                 int                        columns,
                                                 double                     tolerance);
Standard_EXPORT EvolutionData fillBoundary(OcctKernel& kernel, uint32_t wireId, double tolerance);

//! What `patch` built and how it measures: the construction (0 a plane, 1 the Coons patch of
//! four sides, 2 four sides fitted, 3 a trimmed fitted surface, 4 OpenCascade's plate), the
//! degrees and the pole counts in u and v; then per boundary edge the largest gap, angle in
//! degrees, share and difference of curvature; then the largest gap per guide and per point.
struct PatchResult
{
  EvolutionData       evolution;
  std::vector<double> values;
};

//! Controls on this fill; -1 keeps Patch's established automatic construction.
struct PatchControls
{
  int  quad    = -1;
  int  flow    = 0; // Natural, Normal, Next, Adjacent
  int  quality = 0; // Auto, High, Max: minimum fitted spans 6, 10, 16
  bool trim    = true;
  bool open    = false; // Edges that do not close are chained, the gaps closed straight in UV.
};

//! One face inside the closed loop of `edges`. `supports[i]` is the face that edge i is to be
//! smooth with to `orders[i]` (0, 1 or 2 for G0 to G2), or 0 for none; `guides` are edges and
//! `points` coordinates in threes that the face is to pass through. Fails where what it built
//! is measured to miss `tolerance`, `g1Degrees` or `g2Share`. `method` 0 chooses, 1 never
//! takes the plane's shortcut, 2 is OpenCascade's plate.
Standard_EXPORT PatchResult patch(OcctKernel&                  kernel,
                                  const std::vector<uint32_t>& edges,
                                  const std::vector<uint32_t>& supports,
                                  const std::vector<int>&      orders,
                                  const std::vector<uint32_t>& guides,
                                  const std::vector<double>&   points,
                                  double                       tolerance,
                                  double                       g1Degrees,
                                  double                       g2Share,
                                  int                          method,
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
//! placed by the guides. `guideMethod` matters for two guides: 0 turns and scales the profile so
//! that its point on the second guide stays on it (Point), 1 only turns it (Chord), and 2 turns it
//! so that the profile passes through the second guide (Curve). See BRepSkinExt_SweepFramed.cxx.
Standard_EXPORT EvolutionData sweepFramed(OcctKernel&                  kernel,
                                          uint32_t                     profile,
                                          uint32_t                     path,
                                          double                       twist,
                                          double                       scale,
                                          const std::vector<uint32_t>& guides,
                                          int                          guideMethod,
                                          double                       tolerance,
                                          bool                         solid);

//! A skin through two or more profiles along two or more guides, each an edge or a wire, as a
//! Gordon surface; with `solid`, closed flat end profiles are capped. See BRepSkinExt_Gordon.cxx.
Standard_EXPORT EvolutionData skinThrough(OcctKernel&                  kernel,
                                          const std::vector<uint32_t>& sections,
                                          uint32_t                     startVertex,
                                          uint32_t                     endVertex,
                                          bool                         solid,
                                          bool                         ruled);
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
                                            int                          method,
                                            Vec3                         direction,
                                            bool                         both,
                                            bool                         hideOccluded,
                                            int                          completion,
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
//! part only is refused, naming it. `join` as OpenCascade's: 0 arcs, 2 intersections.
Standard_EXPORT EvolutionData offsetSolid(OcctKernel&             kernel,
                                          uint32_t                shape,
                                          double                  distance,
                                          double                  tolerance,
                                          int                     join,
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

Standard_EXPORT EvolutionData sweepSection(OcctKernel&,
                                           uint32_t profile,
                                           uint32_t path,
                                           double   scale,
                                           int      corner,
                                           bool     parallel,
                                           bool     solid,
                                           double   tolerance);
Standard_EXPORT EvolutionData splitFaceEdges(OcctKernel&,
                                             uint32_t                     shape,
                                             uint32_t                     face,
                                             const std::vector<uint32_t>& curves,
                                             double                       tolerance,
                                             const std::vector<int>&      hashes,
                                             int                          bound);

} // namespace BRepExt
