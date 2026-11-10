// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepOffsetExt_Reconstruct.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepPrim_Cylinder.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_FaceInfo.hxx>
#include <BOPDS_PaveBlock.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>

namespace BRepExt
{
namespace
{
gp_Pln outwardPlane(const TopoDS_Face& f)
{
  BRepAdaptor_Surface surface(f);
  gp_Pnt              at;
  gp_Vec              du, dv;
  surface.D1((surface.FirstUParameter() + surface.LastUParameter()) / 2,
             (surface.FirstVParameter() + surface.LastVParameter()) / 2,
             at,
             du,
             dv);
  gp_Dir normal(du.Crossed(dv));
  if (f.Orientation() == TopAbs_REVERSED)
    normal.Reverse();
  return gp_Pln(at, normal);
}

TopoDS_Shape simplify(const TopoDS_Shape& shape, Lineage& lineage)
{
  ShapeUpgrade_UnifySameDomain unify(shape, true, true, false);
  unify.Build();
  lineage.throughHistory(unify.History());
  return unify.Shape();
}

bool usable(const TopoDS_Shape& shape)
{
  if (shape.IsNull() || !TopExp_Explorer(shape, TopAbs_SOLID).More()
      || !BRepCheck_Analyzer(shape).IsValid())
    return false;
  return !intersectsItself(shape);
}
} // namespace

void recordSupportContacts(BRepAlgoAPI_Common& stage,
                           const std::vector<TopoDS_Shape>& sources,
                           const Lineage& facets, const Handle(BRepTools_History)& contacts)
{
  contacts->Merge(stage.History());
  const auto& ds = stage.DSFiller()->DS();
  const auto& builder = *stage.Builder();
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> result;
  TopExp::MapShapes(stage.Shape(), result);
  for (const auto& source : sources)
    for (const auto& facet : facets.of(source))
    {
      const int index = ds.Index(facet);
      if (index < 0 || !ds.HasFaceInfo(index)) continue;
      const auto& info = ds.FaceInfo(index);
      NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> carriers;
      const auto append = [&](int i) {
        if (i < 0) return;
        i = ds.GetSameDomainIndex(i);
        const auto& cell = ds.Shape(i);
        const auto* images = builder.Images().Seek(cell);
        if (images)
          for (const auto& image : *images) carriers.Add(image);
        else carriers.Add(cell);
      };
      for (const auto* blocks : {&info.PaveBlocksIn(), &info.PaveBlocksOn(), &info.PaveBlocksSc()})
        for (const auto& block : *blocks)
        {
          const auto real = ds.RealPaveBlock(block);
          if (real->HasEdge()) append(real->Edge());
          append(real->Pave1().Index());
          append(real->Pave2().Index());
        }
      for (const auto* vertices : {&info.VerticesIn(), &info.VerticesOn(), &info.VerticesSc()})
        for (const int vertex : *vertices) append(vertex);
      for (const auto& carrier : carriers)
      {
        const auto* same = builder.ShapesSD().Seek(carrier);
        const auto& cell = same ? *same : carrier;
        const int found = result.FindIndex(cell);
        if (found && (cell.ShapeType() == TopAbs_EDGE || cell.ShapeType() == TopAbs_VERTEX))
          contacts->AddGenerated(source, result(found));
      }
    }
}

void supportIntersections(const TopoDS_Shape& input, const TopoDS_Shape& result,
                          Lineage& supports, const Handle(BRepTools_History)& contacts)
{
  // The generators here are the constructor's halfspaces, not fitted surfaces.
  // In a convex polyhedron their distinct facets meet in at most one ridge or
  // corner. More than three incident supports allow a non-simple corner to split.
  for (const auto type : {TopAbs_EDGE, TopAbs_VERTEX})
  {
    Adjacency incident;
    TopExp::MapShapesAndUniqueAncestors(input, type, TopAbs_FACE, incident);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> products;
    TopExp::MapShapes(result, type, products);
    for (int i = 1; i <= incident.Extent(); ++i)
    {
      std::vector<TopoDS_Shape> images;
      for (const auto& product : products)
      {
        std::vector<gp_XYZ> normals;
        for (const auto& source : incident(i))
        {
          bool witnessed = false;
          for (const auto* cells : {&contacts->Modified(source), &contacts->Generated(source)})
            for (const auto& carrier : *cells)
            {
              if (carrier.IsSame(product)) witnessed = true;
              for (TopExp_Explorer cell(carrier, type); cell.More(); cell.Next())
                witnessed = witnessed || cell.Current().IsSame(product);
            }
          if (witnessed)
            normals.push_back(BRepAdaptor_Surface(TopoDS::Face(source)).Plane().Axis().Direction().XYZ());
        }
        // Enumerate every independent basis, including non-simple corners. Rank
        // concerns the supplied plane coefficients; no angular cutoff is used.
        bool independent = false;
        for (size_t a = 0; a < normals.size(); ++a)
          for (size_t b = a+1; b < normals.size(); ++b)
          {
            const auto cross = normals[a].Crossed(normals[b]);
            if (type == TopAbs_EDGE) independent = independent || cross.SquareModulus() != 0;
            else
              for (size_t c = b+1; c < normals.size(); ++c)
                independent = independent || cross.Dot(normals[c]) != 0;
          }
        if (independent) images.push_back(product);
      }
      supports.set(incident.FindKey(i), std::move(images));
    }
  }
}

TopoDS_Shape reconstructOffset(const TopoDS_Shape&              solid,
                               const std::vector<TopoDS_Shape>& selected,
                               double                           distance,
                               double                           tolerance,
                               Lineage&                         faces)
{
  if (selected.size() != 1)
    return {};
  const auto          moved = TopoDS::Face(selected[0]);
  BRepAdaptor_Surface surface(moved);
  Adjacency           adjacency;
  TopExp::MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, adjacency);
  TopoDS_Shape result;
  // A full cylindrical cavity bounded only by perpendicular planes may break
  // through surrounding walls. Recut its exact finite analytic cylinder.
  if (surface.GetType() == GeomAbs_Cylinder && moved.Orientation() == TopAbs_REVERSED
      && distance < 0
      && std::abs(surface.LastUParameter() - surface.FirstUParameter() - 2 * M_PI)
           < BRepExt::tolerance::angularResolution)
  {
    const auto c             = surface.Cylinder();
    bool       perpendicular = true;
    for (TopExp_Explorer e(moved, TopAbs_EDGE); e.More(); e.Next())
      for (const auto& f : adjacency.FindFromKey(e.Current()))
        if (!f.IsSame(moved))
        {
          BRepAdaptor_Surface s(TopoDS::Face(f));
          if (s.GetType() != GeomAbs_Plane
              || std::abs(s.Plane().Axis().Direction().Dot(c.Axis().Direction()))
                   < 1 - BRepExt::tolerance::angle)
            perpendicular = false;
        }
    if (perpendicular)
    {
      const auto start =
        c.Location().Translated(gp_Vec(c.Axis().Direction()) * surface.FirstVParameter());
      BRepPrimAPI_MakeCylinder cylinder(gp_Ax2(start, c.Axis().Direction()),
                                       c.Radius() - distance,
                                       surface.LastVParameter() - surface.FirstVParameter());
      const auto tool = cylinder.Shape();
      auto& primitive = cylinder.Cylinder();
      Lineage recut(solid);
      recut.set(moved, {primitive.LateralFace()});
      // The finite cylinder constructor supplies its two circular rims and its
      // periodic seam. Their roles are read in the input parameter domain.
      const double middle = (surface.FirstVParameter() + surface.LastVParameter()) / 2;
      for (TopExp_Explorer it(moved, TopAbs_EDGE); it.More(); it.Next())
      {
        const auto edge = TopoDS::Edge(it.Current());
        BRepAdaptor_Curve2d curve(edge, moved);
        if (BRep_Tool::IsClosed(edge, moved))
        {
          recut.set(edge, {primitive.StartEdge()});
          for (TopExp_Explorer vt(edge, TopAbs_VERTEX); vt.More(); vt.Next())
          {
            const auto vertex = TopoDS::Vertex(vt.Current());
            const bool bottom = BRep_Tool::Parameters(vertex, moved).Y() < middle;
            recut.set(vertex, {bottom ? primitive.BottomStartVertex() : primitive.TopStartVertex()});
          }
        }
        else
        {
          const bool bottom = curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2).Y() < middle;
          recut.set(edge, {bottom ? primitive.BottomEdge() : primitive.TopEdge()});
        }
      }
      Untouched<BRepAlgoAPI_Cut> cut(solid, tool);
      if (cut.IsDone())
      {
        recut.through(cut);
        result = simplify(cut.Shape(), recut);
        faces = std::move(recut);
      }
    }
  }
  if (surface.GetType() == GeomAbs_Plane)
  {
    const auto          plane  = outwardPlane(moved);
    bool                convex = true;
    std::vector<gp_Pln> planes;
    std::vector<TopoDS_Shape> sources;
    for (TopExp_Explorer it(solid, TopAbs_FACE); it.More() && convex; it.Next())
    {
      const auto f = TopoDS::Face(it.Current());
      if (BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane)
      {
        convex = false;
        break;
      }
      auto p = outwardPlane(f);
      for (TopExp_Explorer v(solid, TopAbs_VERTEX); v.More(); v.Next())
        if (gp_Vec(p.Location(), BRep_Tool::Pnt(TopoDS::Vertex(v.Current())))
              .Dot(gp_Vec(p.Axis().Direction()))
            > tolerance * 10)
          convex = false;
      if (f.IsSame(moved))
        p.Translate(gp_Vec(p.Axis().Direction()) * distance);
      planes.push_back(p);
      sources.push_back(f);
    }
    if (convex)
    {
      Bnd_Box b;
      BRepBndLib::Add(solid, b);
      double x0, y0, z0, x1, y1, z1;
      b.Get(x0, y0, z0, x1, y1, z1);
      const double margin = std::abs(distance) + std::sqrt(b.SquareExtent()) + 1;
      result              = BRepPrimAPI_MakeBox(gp_Pnt(x0 - margin, y0 - margin, z0 - margin),
                                   gp_Pnt(x1 + margin, y1 + margin, z1 + margin))
                 .Shape();
      Handle(BRepTools_History) contacts = new BRepTools_History;
      for (size_t i = 0; i < planes.size(); ++i)
      {
        const auto& p = planes[i];
        auto boundary = BRepBuilderAPI_MakeFace(p).Face();
        auto half =
          BRepPrimAPI_MakeHalfSpace(boundary,
                                    p.Location().Translated(-gp_Vec(p.Axis().Direction())))
            .Solid();
        BRepAlgoAPI_Common common(result, half);
        if (!common.IsDone())
          return {};
        faces.set(sources[i], {boundary});
        contacts->AddModified(sources[i], boundary);
        recordSupportContacts(common, sources, faces, contacts);
        faces.through(common);
        result = common.Shape();
      }
      supportIntersections(solid, result, faces, contacts);
      result = simplify(result, faces);
    }
    else if (distance > 0)
    {
      // Extruding a planar face is exactly a face offset only when all
      // its adjacent supports are perpendicular planes. This allows a
      // concave step to disappear; it is not used across angled supports.
      bool perpendicular = true;
      for (TopExp_Explorer e(moved, TopAbs_EDGE); e.More(); e.Next())
        for (const auto& f : adjacency.FindFromKey(e.Current()))
          if (!f.IsSame(moved))
          {
            BRepAdaptor_Surface s(TopoDS::Face(f));
            if (s.GetType() != GeomAbs_Plane
                || std::abs(s.Plane().Axis().Direction().Dot(plane.Axis().Direction()))
                     > BRepExt::tolerance::angle)
              perpendicular = false;
          }
      if (perpendicular)
      {
        BRepPrimAPI_MakePrism prism(moved, gp_Vec(plane.Axis().Direction()) * distance);
        const auto tool = prism.Shape();
        Untouched<BRepAlgoAPI_Fuse> fuse(solid, tool);
        if (fuse.IsDone())
        {
          faces.through(fuse);
          Lineage cap(tool);
          cap.through(fuse);
          ShapeUpgrade_UnifySameDomain unify(fuse.Shape(), true, true, false);
          unify.Build();
          result = unify.Shape();
          faces.throughHistory(unify.History());
          cap.throughHistory(unify.History());
          cap.within(result);
          NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> sources;
          for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
            TopExp::MapShapes(moved, type, sources);
          for (const auto& source : sources)
          {
            auto images = source.ShapeType() == TopAbs_FACE ? std::vector<TopoDS_Shape>{}
                                                          : faces.of(source);
            const auto end = prism.LastShape(source);
            if (!end.IsNull())
            {
              const auto& products = cap.of(end);
              images.insert(images.end(), products.begin(), products.end());
            }
            faces.set(source, std::move(images));
          }
        }
      }
    }
  }
  if (!usable(result))
    return {};
  faces.within(result);
  return solidIfSingle(result);
}
} // namespace BRepExt
