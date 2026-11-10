// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef BRepBlendExt_BoundaryGraph_HeaderFile
#define BRepBlendExt_BoundaryGraph_HeaderFile

#include <TopoDS_Edge.hxx>
#include <map>
#include <optional>
#include <vector>

namespace BRepExt::blending::detail
{
//! Directed occurrences of retained boundary edges. Events identify ports in the
//! sheet parameter domain, not spatial vertices: two ports of a periodic seam
//! may share one spatial vertex. Contact and stop labels do not determine
//! incidence. A contact can occur repeatedly and a component can have any
//! number of edges. Exposure and loop nesting are supplied by the construction.
struct BoundaryGraph
{
  struct Arc
  {
    TopoDS_Edge edge;
    size_t      from, to;
  };

  std::vector<Arc> arcs;

  //! Retained manifold loops, starting at their first stored occurrence. Edge
  //! objects, orientations and occurrence order within a cycle are preserved.
  //! Ambiguous or open incidence is unresolved, never closed by a proximity
  //! threshold or by a choice between outgoing edges.
  std::optional<std::vector<std::vector<TopoDS_Edge>>> loops() const
  {
    std::map<size_t, size_t> outgoing, incoming;
    for (size_t i = 0; i < arcs.size(); ++i)
      if (arcs[i].edge.IsNull() || !outgoing.emplace(arcs[i].from, i).second
          || !incoming.emplace(arcs[i].to, i).second)
        return {};
    if (outgoing.size() != incoming.size())
      return {};
    for (const auto& [port, edge] : outgoing)
      if (incoming.find(port) == incoming.end())
        return {};
    std::vector<std::vector<TopoDS_Edge>> result;
    std::vector<bool>                     used(arcs.size(), false);
    for (size_t i = 0; i < arcs.size(); ++i)
    {
      if (used[i])
        continue;
      std::vector<TopoDS_Edge> loop;
      size_t                   j = i;
      do
      {
        used[j] = true;
        loop.push_back(arcs[j].edge);
        j = outgoing.at(arcs[j].to);
      } while (j != i);
      result.push_back(std::move(loop));
    }
    return result;
  }
};
} // namespace BRepExt::blending::detail
#endif
