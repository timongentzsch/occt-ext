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

#pragma once

// Rolling-ball blends: one blend pipeline built from the exact definition of the blend.
// A blend is defined by its two supports, the offset distances off each, and a section law;
// the ball's centre (the spine) and its feet (the spring curves) are marched exactly, the
// blend surface is approximated from the exact sections with knots at the data's breaks, and
// the blend is attached by cutting the supports back to its contacts and the stopping faces
// back to its ends.

#include <Standard_Macro.hxx>
#include <BRepTools_History.hxx>
#include <TopoDS_Shape.hxx>

#include <utility>
#include <vector>

namespace BRepExt::blending
{

//! Section law: rational conic, curvature-continuous spline, bevel, chord-controlled or
//! three-support round.
enum class Law
{
  Conic   = 0,
  G2      = 1,
  Chamfer = 2,
  Chordal = 3,
  Full    = 4
};

//! Deferred work follows the complete bounded construction pipeline. Hosts with other
//! constructions may request Bounded first and defer further work until those also refuse.
enum class FitWork { Deferred, Bounded, Extended };

struct Request
{
  FitWork fitWork = FitWork::Deferred;
  //! Final contact-arrangement construction, after both existing fit pipelines refuse.
  bool refineContacts = false;
  //! Last construction: stations on the intersection of constant-distance offset supports.
  bool intersectSpines = false;
  const char*               operation = "blendEdges";
  TopoDS_Shape              body;
  std::vector<TopoDS_Shape> edges;
  //! Offset distances of the supports, as a ball's radius is; `second` zero means `first`.
  //! `first` is off `reference` where one is given, else off the first face of each edge.
  double       first = 0, second = 0;
  TopoDS_Shape reference;
  Law          law = Law::Conic;
  //! A chamfer's angle off `reference`, in (0, pi / 2): the setback on the other face is the one
  //! whose bevel makes this angle with the reference face, by the faces' opening at each place
  //! along the edge (`second` is not read). Zero: two distances.
  double angle = 0;
  double rho   = 0.5;
  //! Where not empty, the only edges that tangent ones may be added from (the edges asked are
  //! kept whatever it holds).
  std::vector<TopoDS_Shape> mayAdd;
  //! A size that changes along an edge: for edges[i], where radii[i] is not empty, (place, size)
  //! pairs, the place a fraction of the edge's parameter range from its first parameter to its
  //! last, the size that of `first` there (`second` keeps its ratio to it), interpolated by the
  //! cubic spline with no curvature at its ends (periodic round an edge that closes), as
  //! The reference modeller's variable blends are. Each section is the one of that size in the ball's plane
  //! across the edge. For a round, a conic or a G2 blend of equal distances, or a chamfer; such
  //! edges add no tangent edges.
  std::vector<std::vector<std::pair<double, double>>> radii;
  //! For a full round, the radius it must have all along: the two rounds of this radius whose
  //! contacts meet on the face between, made as two faces meeting where they touch it.
  double meet = 0;
  //! How many replays (replayOverrun) this request is made inside: bounded, so that rounds taken
  //! off and put back that fail the same way end.
  int replays = 0;
};

//! The blended body, made on a copy of the input without its meshes (it shares nothing with the
//! input); what became of every face, edge and vertex of the input (for the face history and the
//! journal); the faces and edges the blends left as they were, each with its copy; and what was
//! made: the pieces, the edges added as tangent to the picked ones and the pieces on exact sheets.
//! Construction intent, independent of the mutable BRep continuity flags. A sharp
//! join remains in the inventory with order zero; it is not a missing contact.
enum class JoinRole
{
  Contact, UsedSupport, Joint, SupportBreak, Fold, Roll, Corner, Miter, Stop,
  Trim, Ridge, Split
};

struct Join
{
  TopoDS_Shape edge;
  JoinRole role;
  int order; //!< Required geometric continuity: 0, 1 or 2; not a measured result.
};

struct Blended
{
  TopoDS_Shape                                       shape;
  Handle(BRepTools_History)                          history = new BRepTools_History;
  std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> unchanged;
  std::vector<Join> joins;
};

Standard_EXPORT Blended blend(const Request& request);

//! The kernel's stored body `shape` blended at its `edges` and checked valid: what the blend's
//! two entry points, with a journal session and without, share. `law` is a Law's value and
//! `reference` 0 for none.
//! A result the blend may return for this input: valid, closed where the input was, and every edge
//! that claims it curvature continuous to within limit::seamBend of a blend of `size`. Refused by
//! name otherwise.
Standard_EXPORT void accepted(const char*         operation,
                              const TopoDS_Shape& input,
                              const TopoDS_Shape& result,
                              double              size,
                              const Blended*      made = nullptr);

} // namespace BRepExt::blending
