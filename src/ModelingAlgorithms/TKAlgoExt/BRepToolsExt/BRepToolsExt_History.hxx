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

// Face history in upstream's EvolutionData format, so TypeScript reads our results exactly like
// the facade's own *WithHistory results:
//   modified / generated: flat records [inputFaceHash, count, resultFaceHash...]
//   deleted:              [inputFaceHash...]
// Hashes are TopTools_ShapeMapHasher % hashUpperBound, the same value OcctKernel::hashCode gives.
// Only input faces whose hash the caller listed are reported.

#include <BRepToolsExt_Checks.hxx>

#include <BRepTools_History.hxx>
#include <NCollection_List.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include <unordered_set>

namespace BRepExt
{

class FaceHistory
{
public:
  FaceHistory(const char* operation, const std::vector<int>& inputFaceHashes, int hashUpperBound)
      : tracked_(inputFaceHashes.begin(), inputFaceHashes.end()),
        bound_(hashUpperBound)
  {
    // Every hash is a remainder by this bound, and a zero divisor traps.
    if (hashUpperBound <= 0)
      fail(operation, invalid_input, "hashUpperBound must be positive");
  }

  int hash(const TopoDS_Shape& shape) const
  {
    return static_cast<int>(TopTools_ShapeMapHasher{}(shape) % static_cast<size_t>(bound_));
  }

  bool tracks(const TopoDS_Shape& face) const { return tracked_.count(hash(face)) != 0; }

  void modified(const TopoDS_Shape& face, const NCollection_List<TopoDS_Shape>& successors)
  {
    record(evolution_.modified, face, successors);
  }

  void generated(const TopoDS_Shape& face, const NCollection_List<TopoDS_Shape>& successors)
  {
    record(evolution_.generated, face, successors);
  }

  void deleted(const TopoDS_Shape& face) { evolution_.deleted.push_back(hash(face)); }

  //! Record what an OCCT maker or BRepTools_History reports for
  //! each tracked face of `input`. Same reading as upstream's appendEvolution.
  template <typename Maker>
  void collect(Maker& maker, const TopoDS_Shape& input)
  {
    for (TopExp_Explorer faces(input, TopAbs_FACE); faces.More(); faces.Next())
    {
      const TopoDS_Shape& face = faces.Current();
      if (!tracks(face))
        continue;
      modified(face, maker.Modified(face));
      generated(face, maker.Generated(face));
      if (deletedBy(maker, face))
        deleted(face);
    }
  }

  // Compose two topology-changing stages; never publish intermediate hashes.
  template <typename First, typename Second>
  void collectComposed(First& first, Second& second, const TopoDS_Shape& input)
  {
    collectComposed(first, second, input, [](const TopoDS_Shape&, NCollection_List<TopoDS_Shape>&) {
    });
  }

  //! The same, with `extend(face, successors)` adding to an input face's successors what the
  //! stages do not report themselves.
  template <typename First, typename Second, typename Extend>
  void collectComposed(First& first, Second& second, const TopoDS_Shape& input, Extend extend)
  {
    for (TopExp_Explorer it(input, TopAbs_FACE); it.More(); it.Next())
    {
      const auto original = it.Current();
      if (!tracks(original))
        continue;
      auto middle = first.Modified(original);
      if (middle.IsEmpty() && !deletedBy(first, original))
        middle.Append(original);
      NCollection_List<TopoDS_Shape> final;
      for (const auto& face : middle)
      {
        const auto& changed = second.Modified(face);
        if (changed.IsEmpty() && !deletedBy(second, face))
          final.Append(face);
        else
          for (const auto& successor : changed)
            final.Append(successor);
      }
      extend(original, final);
      if (final.IsEmpty())
        deleted(original);
      else
        modified(original, final);
    }
  }

  EvolutionData finish(uint32_t resultId)
  {
    evolution_.resultId = resultId;
    return std::move(evolution_);
  }

private:
  template <typename Maker>
  static bool deletedBy(Maker& maker, const TopoDS_Shape& face)
  {
    return maker.IsDeleted(face);
  }

  static bool deletedBy(BRepTools_History& history, const TopoDS_Shape& face)
  {
    return history.IsRemoved(face);
  }

  void record(std::vector<int>&                     into,
              const TopoDS_Shape&                   face,
              const NCollection_List<TopoDS_Shape>& successors)
  {
    // OCCT also lists the edges and vertices a face gave rise to; only faces are history here.
    const size_t header = into.size();
    into.push_back(hash(face));
    into.push_back(0);
    for (const TopoDS_Shape& successor : successors)
    {
      if (successor.ShapeType() != TopAbs_FACE)
        continue;
      into.push_back(hash(successor));
      into[header + 1]++;
    }
    if (into[header + 1] == 0)
      into.resize(header);
  }

  std::unordered_set<int> tracked_;
  int                     bound_;
  EvolutionData           evolution_;
};

} // namespace BRepExt
