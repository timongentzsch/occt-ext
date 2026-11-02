// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The result check every operation makes (BRepToolsExt_Checks.hxx): it must agree with
// OpenCascade's full checker on what a body is, and its incremental form must look at what the
// operation made, as documented, and not at what it was given.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>

#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

// A "solid" built on a shell of five of the six faces of a box: closed nowhere, so any checker must
// reject it.
TopoDS_Solid openSolid()
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(5, 4, 3).Shape();
  BRep_Builder       builder;
  TopoDS_Shell       shell;
  builder.MakeShell(shell);
  int count = 0;
  for (TopExp_Explorer it(box, TopAbs_FACE); it.More() && count < 5; it.Next(), ++count)
    builder.Add(shell, it.Current());
  return BRepBuilderAPI_MakeSolid(shell).Solid();
}

} // namespace

// Valid and invalid bodies get the verdict of the full checker, and the report is empty exactly
// when the body is valid.
TEST(BRepToolsExtChecks, AgreesWithTheFullCheckerOnValidAndInvalidBodies)
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(5, 4, 3).Shape();
  EXPECT_TRUE(BRepCheck_Analyzer(box).IsValid());
  EXPECT_TRUE(BRepExt::invalidity(box).empty());
  const TopoDS_Solid open = openSolid();
  EXPECT_FALSE(BRepCheck_Analyzer(open).IsValid());
  EXPECT_FALSE(BRepExt::invalidity(open).empty());
}

// For the result of an operation, the incremental check finds what is wrong in what the operation
// made (a shell that is not closed is new, so it is reported) and agrees with the full check on a
// valid result of a real fillet, whose faces are partly new and partly shared with the input.
TEST(BRepToolsExtChecks, IncrementalCheckLooksAtWhatTheOperationMade)
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(5, 4, 3).Shape();
  EXPECT_FALSE(BRepExt::invalidity(openSolid(), box).empty());
  BRepFilletAPI_MakeFillet maker(box);
  for (TopExp_Explorer it(box, TopAbs_EDGE); it.More(); it.Next())
    maker.Add(0.5, TopoDS::Edge(it.Current()));
  maker.Build();
  ASSERT_TRUE(maker.IsDone());
  const TopoDS_Shape filleted = maker.Shape();
  EXPECT_TRUE(BRepCheck_Analyzer(filleted).IsValid());
  EXPECT_TRUE(BRepExt::invalidity(filleted, box).empty());
}
