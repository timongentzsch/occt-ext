// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#include <Standard.hxx>
#include <Message_ProgressRange.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <gp_Vec.hxx>
#include <functional>
#include <map>
#include <vector>
namespace GeomAlgoExt {
struct BlendContactJet { gp_Vec du, dv, duu, dvv, duv, n; };
struct BlendFitResult {
 Handle(Geom_BSplineSurface) blend;
 Handle(Geom2d_BSplineCurve) contact[2];
};
//! Exact section source and caller-owned law/tolerance policy. No geometric law defaults.
struct BlendFitInput {
 int degree, sectionDegree, sectionCount;
 std::vector<double> sectionKnots;
 std::vector<int> sectionMultiplicities;
 std::map<double,int> knots;
 double scale, tolerance, knotMerge, degenerate, bendTolerance;
 int bendRefinements;
 bool curvature, contactPilot, keepFirst, keepLast, cliff[2];
 std::function<std::vector<double>(double)> row;
 std::function<std::vector<double>(double,int)> terminalRow;
 std::function<std::map<double,int>()> discoveredBreaks;
 std::function<int(double)> regime;
 std::function<void(std::map<double,int>&,double,int)> addBreak;
 std::function<bool(int,double,double,BlendContactJet&)> contactJet;
 std::function<gp_Pnt(int,double,double)> contactPoint;
};
//! A single original refinement iteration. The input must outlive this state.
//! This numerical step has no internal cancellation points; callers own the outer loop.
class BlendSectionFitter {
 const BlendFitInput& myInput;
 std::map<double,int> myKnots;
 std::vector<double> myWorst;
 int myBendRounds=0;
 bool myStopped=false;
 BlendFitResult myResult;
public:
 Standard_EXPORT BlendSectionFitter(const BlendFitInput&,std::map<double,int> knots);
 Standard_EXPORT bool Advance(int refinement);
 bool IsStopped() const {return myStopped;}
 const BlendFitResult& Result() const {return myResult;}
};
Standard_EXPORT BlendFitResult fitBlendSections(const BlendFitInput&, const Message_ProgressRange&);
}
