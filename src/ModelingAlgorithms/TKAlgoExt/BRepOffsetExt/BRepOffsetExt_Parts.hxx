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

// What shells, solid offsets and face moves share: which faces an offset cannot keep, what
// became of each face through several stages, and the measurements a result is held to.

#include <Standard_Macro.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_History.hxx>

#include <BRepTools_History.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopoDS_Face.hxx>

#include <vector>
#include <limits>
#include <cmath>

namespace BRepExt
{

// --- Tolerances of this subject, each with where its value is from -------------------------

//! A face counts as gone at a point where the offset reaches this share of its radius of
//! curvature. OpenCascade's blends on curved supports are approximations to 1e-4 of their
//! radius (ChFi3d's `tolapp3d`), so a wall equal to such a radius has to count as reaching it.
inline constexpr double collapseSlack = 1e-4;

//! How many times the caller's tolerance a wall may be thinner than asked where it is
//! measured. The offset of a spline face is itself fitted to that tolerance, and the distance
//! to the input is found by an iteration that stops at 1e-7.
inline constexpr double clearanceSlack = 100.0;

//! By what share of the wall the tool that opens a shell stands out of the face it opens. It
//! only has to be far above the Boolean's own fuzz (1e-7) and small against the part, so that
//! it reaches nothing that overhangs the opening.
inline constexpr double openingOvershoot = 0.05;

// --- Faces that an offset cannot keep ------------------------------------------------------

//! The faces of a body that an offset by `distance` along the outward normal (negative:
//! inwards) cannot keep: `vanish` where the whole face is bent tighter than the distance on
//! the side it moves to, `partial` where only a part of it is.
struct OffsetFate
{
  std::vector<TopoDS_Shape> vanish, partial;
};

Standard_EXPORT OffsetFate classifyForOffset(const TopoDS_Shape& body, double distance);

//! Whether faces of `shape` cross each other, which BRepCheck_Analyzer does not look for.
Standard_EXPORT bool intersectsItself(const TopoDS_Shape& shape);

//! Points inside `face`, away from its boundary: up to `perSide` squared of them.
Standard_EXPORT std::vector<gp_Pnt> pointsInside(const TopoDS_Face& face, int perSide);

//! Coverage and the least measured distance are separate: an empty or unresolved
//! region cannot masquerade as infinite clearance.
struct Clearance
{
  double minimum = std::numeric_limits<double>::infinity();
  size_t covered = 0, expected = 0;
  bool resolved = true;
  std::string issue;
  bool satisfies(double required) const
  {
    return resolved && expected > 0 && covered == expected
           && std::isfinite(minimum) && minimum >= required;
  }
};

//! Whole trimmed-face distances. Near openings, subdivide the trimmed domain
//! until each region is measured clear or proved inside the opening exclusion.
//! A work limit leaves coverage unresolved; it never counts as a measurement.
Standard_EXPORT Clearance clearance(const std::vector<TopoDS_Shape>& faces,
                                    const TopoDS_Shape&              body,
                                    const std::vector<TopoDS_Shape>& opened = {},
                                    double wall = 0.0,
                                    double required = 0.0);

//! Push-pull by how much the far side of a face may lean from the direction it moves in and
//! still count as running along it, as a cosine: 1e-7 is an angle of 1e-7 radians, which
//! over a part of a metre is the 1e-7 that OpenCascade takes two points for one at.
inline constexpr double alongSlack = 1e-7;

//! Whether `a` and `b` are tangent along their common `edge`.
Standard_EXPORT bool tangentAlong(const TopoDS_Edge& edge,
                                  const TopoDS_Face& a,
                                  const TopoDS_Face& b);

// --- What became of each face --------------------------------------------------------------

//! The faces that each face of an input has become, carried through any number of stages.
class Lineage
{
public:
  explicit Lineage(const TopoDS_Shape& input);

  //! A stage that reports through Modified and IsDeleted: a Boolean, a defeaturing.
  template <class Stage>
  void through(Stage& stage)
  {
    step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
      const NCollection_List<TopoDS_Shape>& next = stage.Modified(face);
      if (next.IsEmpty() && !stage.IsDeleted(face))
        out.push_back(face);
      for (const TopoDS_Shape& successor : next)
        if (successor.ShapeType() == TopAbs_FACE)
          out.push_back(successor);
    });
  }

  //! A stage whose faces are images of the input's, reported through Generated: an offset.
  template <class Stage>
  void images(Stage& stage)
  {
    step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
      for (const TopoDS_Shape& successor : stage.Generated(face))
        if (successor.ShapeType() == TopAbs_FACE)
          out.push_back(successor);
    });
  }

  //! A stage that files some successors under Modified and others under Generated: a draft.
  template <class Stage>
  void throughBoth(Stage& stage)
  {
    step([&](const TopoDS_Shape& face, std::vector<TopoDS_Shape>& out) {
      const size_t before = out.size();
      for (const TopoDS_Shape& successor : stage.Modified(face))
        if (successor.ShapeType() == TopAbs_FACE)
          out.push_back(successor);
      for (const TopoDS_Shape& successor : stage.Generated(face))
        if (successor.ShapeType() == TopAbs_FACE)
          out.push_back(successor);
      if (out.size() == before && !stage.IsDeleted(face))
        out.push_back(face);
    });
  }

  void throughHistory(const Handle(BRepTools_History)& history);

  //! What `face` has become, said outright: where a stage has no report of its own.
  void set(const TopoDS_Shape& face, std::vector<TopoDS_Shape> now);
  //! Then through a later operation, of which `next` is the lineage: each face becomes
  //! what `next` has it become; one that `next` does not know stays.
  void then(const Lineage& next);

  //! Keep only what is a face of `result`, as `result` holds it.
  void within(const TopoDS_Shape& result);

  const std::vector<TopoDS_Shape>& of(const TopoDS_Shape& face) const;
  //! Every face that descends from any face of the input.
  std::vector<TopoDS_Shape> all() const;

  //! Into the history: as what each tracked face was modified to (a face that is itself is
  //! left out, one with nothing left is deleted), or as what it generated.
  void modified(FaceHistory& history) const;
  void generated(FaceHistory& history) const;

private:
  template <class Next>
  void step(Next&& next)
  {
    for (auto& faces : now_)
    {
      std::vector<TopoDS_Shape> out;
      for (const TopoDS_Shape& face : faces)
        next(face, out);
      faces = std::move(out);
    }
  }

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> input_;
  std::vector<std::vector<TopoDS_Shape>>                        now_;
};

//! `body` with each face on a surface C0 in its parameters divided along those knots, as
//! OpenCascade's offset refuses such a surface however smooth it is (BRepOffset_C0Geometry). The
//! geometry is the same; `body` itself where it has no such face.
Standard_EXPORT TopoDS_Shape dividedAtC0(const TopoDS_Shape& body, Lineage& faces);

} // namespace BRepExt
