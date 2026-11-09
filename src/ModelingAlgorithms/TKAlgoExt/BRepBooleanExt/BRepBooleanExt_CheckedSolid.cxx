// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepToolsExt_History.hxx>
#include <GeomAlgoExt_Timing.hxx>

#include <BRepBuilderAPI_ModifyShape.hxx>
#include <BRepTools_CopyModification.hxx>
#include <Geom_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <unordered_map>
#include <algorithm>
#include <NCollection_Map.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>

// The kernel's Booleans, as the commands call them. All of them are selections of the cells of
// one General Fuse, verified against the operands as they were given: see BRepBooleanExt_Cells.cxx.

namespace BRepExt
{
namespace
{

// Preserve exact surface sharing while constructing the copied topology. The
// modifier attaches copied edge pcurves to these same owned surface handles.
class SharedSurfaceModification : public BRepTools_CopyModification
{
  std::unordered_map<const Geom_Surface*,
                     std::vector<std::pair<TopLoc_Location, Handle(Geom_Surface)>>>
    surfaces;

public:
  SharedSurfaceModification()
      : BRepTools_CopyModification(true, false)
  {
  }

  bool NewSurface(const TopoDS_Face&    face,
                  Handle(Geom_Surface)& surface,
                  TopLoc_Location&      location,
                  double&               tolerance,
                  bool&                 reverseWires,
                  bool&                 reverseFace) override
  {
    const bool      changed = BRepTools_CopyModification::NewSurface(face,
                                                                surface,
                                                                location,
                                                                tolerance,
                                                                reverseWires,
                                                                reverseFace);
    TopLoc_Location originalLocation;
    const auto      original = BRep_Tool::Surface(face, originalLocation);
    if (original.IsNull())
      return changed;
    auto&      group = surfaces[original.get()];
    const auto known = std::find_if(group.begin(), group.end(), [&](const auto& entry) {
      return entry.first.IsEqual(originalLocation);
    });
    if (known == group.end())
      group.emplace_back(originalLocation, surface);
    else
      surface = known->second;
    return changed;
  }
};

class SharedSurfaceCopy : public BRepBuilderAPI_ModifyShape
{
public:
  explicit SharedSurfaceCopy(const TopoDS_Shape& shape)
      : BRepBuilderAPI_ModifyShape(shape, new SharedSurfaceModification)
  {
  }
};

BooleanRequest requestOf(OcctKernel&                  k,
                         const char*                  operation,
                         Select                       select,
                         const std::vector<uint32_t>& sources,
                         const std::vector<uint32_t>& tools,
                         int                          glue,
                         double                       fuzzy)
{
  BooleanRequest request;
  request.operation = operation;
  request.select    = select;
  request.glue      = glue;
  request.fuzzy     = fuzzy;
  for (auto id : sources)
    request.sources.push_back(get(k, id));
  for (auto id : tools)
    request.tools.push_back(get(k, id));
  return request;
}

//! Faces within `angle` of one surface made one, as BRepAlgoAPI's SimplifyResult does, with
//! what each operand's face became carried through. On a copy: the result shares faces with
//! the operands, and the unifying changes what it works on in place (see `reunite`).
} // namespace

void simplifyBooleanResult(BooleanResult&                   result,
                           const BooleanRequest&            request,
                           double                           angle,
                           const std::vector<TopoDS_Shape>& kept)
{
  SharedSurfaceCopy                                      copier(result.shape);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> present;
  for (TopExp_Explorer it(result.shape, TopAbs_FACE); it.More(); it.Next())
    present.Add(it.Current());
  ShapeUpgrade_UnifySameDomain unify(copier.Shape(), true, true, false);
  for (const auto& shape : kept)
    unify.KeepShape(copier.ModifiedShape(shape));
  unify.SetAngularTolerance(angle);
  unify.Build();
  if (!BRepCheck_Analyzer(unify.Shape()).IsValid())
    return; // the result stays as it was, which is checked
  for (const auto* group : {&request.sources, &request.tools})
    for (const auto& operand : *group)
      for (TopExp_Explorer it(operand, TopAbs_FACE); it.More(); it.Next())
      {
        NCollection_List<TopoDS_Shape> before, after;
        if (const auto* known = result.faces.Seek(it.Current()))
          before = *known;
        else
          before.Append(it.Current());
        for (const auto& original : before)
        {
          if (!present.Contains(original))
            continue;
          const TopoDS_Shape face = copier.ModifiedShape(original);
          const auto&        now  = unify.History()->Modified(face);
          if (now.IsEmpty() && !unify.History()->IsRemoved(face))
            after.Append(face);
          for (const auto& successor : now)
            after.Append(successor);
        }
        // Every face of the copy is another than the one it was copied from.
        result.faces.Bind(it.Current(), after);
      }
  // Keep the topology journal in step with the public face map. Callers
  // composing another operation must never receive intermediate copy hashes.
  Handle(BRepTools_History)                                     step = new BRepTools_History;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> original, final;
  TopExp::MapShapes(result.shape, original);
  TopExp::MapShapes(unify.Shape(), final);
  for (int i = 1; i <= original.Extent(); ++i)
  {
    const auto& before = original(i);
    if (!BRepTools_History::IsSupportedType(before))
      continue;
    const auto  copy     = copier.ModifiedShape(before);
    bool        retained = false;
    const auto& changed  = unify.History()->Modified(copy);
    if (changed.IsEmpty() && !unify.History()->IsRemoved(copy) && final.Contains(copy))
    {
      step->AddModified(before, copy);
      retained = true;
    }
    for (const auto& now : changed)
      if (final.Contains(now))
      {
        step->AddModified(before, now);
        retained = true;
      }
    for (const auto& now : unify.History()->Generated(copy))
      if (final.Contains(now))
        step->AddGenerated(before, now);
    if (!retained)
      step->Remove(before);
  }
  result.history->Merge(step);
  result.shape = unify.Shape();
}

namespace
{

EvolutionData finish(OcctKernel&             k,
                     BooleanResult&          result,
                     const BooleanRequest&   request,
                     const std::vector<int>& hashes,
                     int                     bound)
{
  BooleanPhase phase("history");
  FaceHistory  history(request.operation, hashes, bound);
  for (const auto* group : {&request.sources, &request.tools})
    for (const auto& operand : *group)
      history.collect(result, operand);
  return history.finish(store(k, result.shape));
}

} // namespace

EvolutionData booleanChecked(OcctKernel&                  k,
                             int                          operation,
                             const std::vector<uint32_t>& sources,
                             const std::vector<uint32_t>& tools,
                             int                          glue,
                             double                       fuzzy,
                             double                       angle,
                             const std::vector<int>&      hashes,
                             int                          bound)
{
  if (operation < 0 || operation > 2)
    fail("boolean", invalid_input, "unknown operation");
  auto request =
    requestOf(k, "booleanChecked", static_cast<Select>(operation), sources, tools, glue, fuzzy);
  auto result = booleanCells(request);
  if (angle > 0)
    simplifyBooleanResult(result, request, angle);
  return finish(k, result, request, hashes, bound);
}

EvolutionData booleanCells(OcctKernel&                  k,
                           int                          select,
                           const std::vector<uint32_t>& sources,
                           const std::vector<uint32_t>& tools,
                           double                       fuzzy,
                           const std::vector<int>&      hashes,
                           int                          bound)
{
  if (select < 0 || select > 6)
    fail("boolean", invalid_input, "unknown selection of cells");
  auto request =
    requestOf(k, "booleanCells", static_cast<Select>(select), sources, tools, 0, fuzzy);
  auto result = booleanCells(request);
  return finish(k, result, request, hashes, bound);
}

std::vector<Rung> booleanReport(OcctKernel&                  k,
                                int                          select,
                                const std::vector<uint32_t>& sources,
                                const std::vector<uint32_t>& tools,
                                double                       fuzzy)
{
  if (select < 0 || select > 6)
    fail("boolean", invalid_input, "unknown selection of cells");
  auto request = requestOf(k, "boolean", static_cast<Select>(select), sources, tools, 0, fuzzy);
  std::vector<Rung> rungs;
  try
  {
    booleanCells(request, &rungs);
  }
  catch (const Failure&)
  {
    if (rungs.empty())
      throw; // refused before any rung: the request itself is wrong
  }
  return rungs;
}

} // namespace BRepExt
