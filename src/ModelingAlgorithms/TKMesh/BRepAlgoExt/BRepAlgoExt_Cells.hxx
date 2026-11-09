// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Booleans as selections of the cells of one General Fuse.
//
// All operands go through one BOPAlgo_CellsBuilder, which leaves them as they were. What it
// returns is a partition: every operand is cut into cells (solids into solids, sheets into
// faces, curves into edges), and a cell that two operands share exists once. A union, a
// difference, an intersection and a slice are then choices among the same cells, so they
// cannot disagree with each other.
//
// Before anything is chosen the partition is checked against the ORIGINAL operands
// (see `verify` in BRepBooleanExt_Cells.cxx), and what cannot be verified is refused.

#include <Standard_Macro.hxx>
#include <BRepTools_History.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Shape.hxx>

#include <string>
#include <vector>

namespace BRepExt
{

//! Which cells make the result. `sources` are what is worked on, `tools` what it is worked with.
enum class Select
{
  Union        = 0, // every cell
  Difference   = 1, // cells of the sources outside every tool
  Intersection = 2, // cells of the sources inside some tool (OpenCascade's reading of Common)
  Slice        = 3, // every cell of the sources, each for itself
  Imprint      = 4, // the sources whole again, their faces left divided where the tools met them
  Common       = 5, // cells inside every operand
  Region       = 6, // each bounded region enclosed by all operands
  Behind       = 7, // cells of the sources behind the sheet tools: on the side their faces'
                    // normals point away from (oriented supporting surface)
  Front = 8,        // cells of the sources in front of the sheet tools
  Auto  = 9,        // Difference where a tool reaches inside a source, else Union where they
                    // share a face, else nothing: the choice Extrude makes from the clash
};

//! How a Boolean glues coincident operands; the values are BOPAlgo_GlueEnum's.
enum class BooleanGlue
{
  Off   = 0,
  Shift = 1,
  Full  = 2
};

//! What one attempt of the ladder (plain, then fuzzy) did; these make up its timing table.
struct BooleanAttempt
{
  std::string name;         // "plain", "fuzzy"
  double      fuzzy    = 0; // the fuzzy value it ran at
  bool        accepted = false;
  std::string refusal;      // why not, if not
  int         cells    = 0; // cells of the partition
  int         interior = 0; // classifications of a point inside a cell
  int         probes   = 0; // classifications of a point on a cell's boundary
  double      residual =
    0; // what the cells of an operand and the operand differ by, where that is nearest its limit
  double slack     = 0; // the limit: what it was allowed to be
  double tolerance = 0; // largest tolerance of an edge or face of the partition
  double fuseMs = 0, checkMs = 0, buildMs = 0;
  // Of checkMs: finding a point inside each cell, classifying points, integrating faces.
  double insideMs = 0, classifyMs = 0, measureMs = 0;
};

struct BooleanResult
{
  TopoDS_Shape shape;
  //! What Auto chose: Difference, Union, or None when the tools neither reach into nor touch
  //! the sources (then `shape` is null). Other selections are what was asked.
  enum class Decided
  {
    Asked,
    Difference,
    Union,
    None
  } decided                         = Decided::Asked;
  Handle(BRepTools_History) history = new BRepTools_History;
  //! For each face of an operand that is not in the result as it was: what it became. An
  //! empty list is a face that is gone. A face that is not bound is in the result untouched.
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> faces;
  std::vector<BooleanAttempt> attempts;
  //! Exact copies untouched by reunion, for display reuse; never geometry matching.
  std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> displayCopies;

  // What BRepToolsExt_History.hxx's FaceHistory::collect asks of a maker.
  const NCollection_List<TopoDS_Shape>& Modified(const TopoDS_Shape& face) const
  {
    const auto* found = faces.Seek(face);
    return found ? *found : none_;
  }

  const NCollection_List<TopoDS_Shape>& Generated(const TopoDS_Shape&) const { return none_; }

  bool IsDeleted(const TopoDS_Shape& face) const
  {
    const auto* found = faces.Seek(face);
    return found && found->IsEmpty();
  }

private:
  NCollection_List<TopoDS_Shape> none_;
};

struct BooleanRequest
{
  const char*               operation = "boolean"; // the name failures carry
  Select                    select    = Select::Union;
  std::vector<TopoDS_Shape> sources, tools;
  bool                      cancellable = false;
  bool        subShapeHistory = false; // also the history of edges and vertices, not only faces
  BooleanGlue glue            = BooleanGlue::Off;
  double      fuzzy           = 0; // the caller's own fuzzy value for the first rung; 0 for none
};

//! Run the ladder: plain, then at a fuzzy value; the first rung whose partition is verified
//! gives the result. Fails with no_solution, saying why each rung was refused, when none is.
//! `attempts`, if given, receives what every attempt did, also when the call fails.
Standard_EXPORT BooleanResult booleanCells(const BooleanRequest&        request,
                                           std::vector<BooleanAttempt>* attempts = nullptr);

//! Copy and unite same-domain result faces, composing operand face history. Kept result
//! edges/vertices retain intentional boundaries; invalid unification leaves result unchanged.
Standard_EXPORT void simplifyBooleanResult(BooleanResult&                   result,
                                           const BooleanRequest&            request,
                                           double                           angle,
                                           const std::vector<TopoDS_Shape>& kept = {});

} // namespace BRepExt
