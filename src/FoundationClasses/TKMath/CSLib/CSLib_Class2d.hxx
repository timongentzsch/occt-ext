// Created on: 1995-03-08
// Created by: Laurent BUCHARD
// Copyright (c) 1995-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
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

#ifndef _CSLib_Class2d_HeaderFile
#define _CSLib_Class2d_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <gp_Pnt2d.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_DynamicArray.hxx>

class gp_Pnt2d;

//! Low-level algorithm for 2D point-in-polygon classification.
//!
//! This class determines whether a 2D point lies inside, outside, or on the boundary
//! of a closed polygon. It uses a ray-casting algorithm where a horizontal ray
//! from the test point is extended to infinity, and the number of polygon edge
//! crossings determines the classification.
//!
//! The polygon is internally normalized to [0,1] x [0,1] domain for numerical stability.
//!
//! @note This class was moved from package BRepTopAdaptor.
class CSLib_Class2d
{
public:
  //! Classification result for point-in-polygon tests.
  enum Result
  {
    Result_Inside    = 1,  //!< Point is strictly inside the polygon
    Result_Outside   = -1, //!< Point is strictly outside the polygon
    Result_Uncertain = 0   //!< Point is on boundary or classification is uncertain
  };

  DEFINE_STANDARD_ALLOC

  //! Default constructor. Creates an empty classifier.
  CSLib_Class2d() = default;

  //! Constructs a 2D classifier from an array of polygon vertices.
  //!
  //! The polygon is automatically closed (no need to repeat the first point at the end).
  //! Points are normalized internally to the UV bounds for numerical stability.
  //!
  //! @param[in] thePnts2d Array of polygon vertices (minimum 3 points required)
  //! @param[in] theTolU   Tolerance in U direction for boundary detection
  //! @param[in] theTolV   Tolerance in V direction for boundary detection
  //! @param[in] theUMin   Minimum U bound of the polygon domain
  //! @param[in] theVMin   Minimum V bound of the polygon domain
  //! @param[in] theUMax   Maximum U bound of the polygon domain
  //! @param[in] theVMax   Maximum V bound of the polygon domain
  Standard_EXPORT CSLib_Class2d(const NCollection_Array1<gp_Pnt2d>& thePnts2d,
                                double                              theTolU,
                                double                              theTolV,
                                double                              theUMin,
                                double                              theVMin,
                                double                              theUMax,
                                double                              theVMax);

  //! Constructs a 2D classifier from a sequence of polygon vertices.
  //!
  //! Same as the array constructor but accepts a sequence for convenience.
  //!
  //! @param[in] thePnts2d Sequence of polygon vertices (minimum 3 points required)
  //! @param[in] theTolU   Tolerance in U direction for boundary detection
  //! @param[in] theTolV   Tolerance in V direction for boundary detection
  //! @param[in] theUMin   Minimum U bound of the polygon domain
  //! @param[in] theVMin   Minimum V bound of the polygon domain
  //! @param[in] theUMax   Maximum U bound of the polygon domain
  //! @param[in] theVMax   Maximum V bound of the polygon domain
  Standard_EXPORT CSLib_Class2d(const NCollection_Sequence<gp_Pnt2d>& thePnts2d,
                                double                                theTolU,
                                double                                theTolV,
                                double                                theUMin,
                                double                                theVMin,
                                double                                theUMax,
                                double                                theVMax);

  //! Constructs a 2D classifier from a vector of polygon vertices.
  //!
  //! Same as the array constructor but accepts a vector for convenience.
  //!
  //! @param[in] thePnts2d Vector of polygon vertices (minimum 3 points required)
  //! @param[in] theTolU   Tolerance in U direction for boundary detection
  //! @param[in] theTolV   Tolerance in V direction for boundary detection
  //! @param[in] theUMin   Minimum U bound of the polygon domain
  //! @param[in] theVMin   Minimum V bound of the polygon domain
  //! @param[in] theUMax   Maximum U bound of the polygon domain
  //! @param[in] theVMax   Maximum V bound of the polygon domain
  Standard_EXPORT CSLib_Class2d(const NCollection_DynamicArray<gp_Pnt2d>& thePnts2d,
                                double                                    theTolU,
                                double                                    theTolV,
                                double                                    theUMin,
                                double                                    theVMin,
                                double                                    theUMax,
                                double                                    theVMax);

  //! Move constructor.
  CSLib_Class2d(CSLib_Class2d&& theOther) noexcept
      : myPnts2dX(std::move(theOther.myPnts2dX)),
        myPnts2dY(std::move(theOther.myPnts2dY)),
        myBucketStart(std::move(theOther.myBucketStart)),
        myBucketEdges(std::move(theOther.myBucketEdges)),
        myBucketCount(theOther.myBucketCount),
        myTolU(theOther.myTolU),
        myTolV(theOther.myTolV),
        myPointsCount(theOther.myPointsCount),
        myUMin(theOther.myUMin),
        myVMin(theOther.myVMin),
        myUMax(theOther.myUMax),
        myVMax(theOther.myVMax)
  {
  }

  //! Move assignment operator.
  CSLib_Class2d& operator=(CSLib_Class2d&& theOther) noexcept
  {
    if (this != &theOther)
    {
      myPnts2dX     = std::move(theOther.myPnts2dX);
      myPnts2dY     = std::move(theOther.myPnts2dY);
      myBucketStart = std::move(theOther.myBucketStart);
      myBucketEdges = std::move(theOther.myBucketEdges);
      myBucketCount = theOther.myBucketCount;
      myTolU        = theOther.myTolU;
      myTolV        = theOther.myTolV;
      myPointsCount = theOther.myPointsCount;
      myUMin        = theOther.myUMin;
      myVMin        = theOther.myVMin;
      myUMax        = theOther.myUMax;
      myVMax        = theOther.myVMax;
    }
    return *this;
  }

  //! Classifies a point relative to the polygon.
  //!
  //! @param[in] thePoint The 2D point to classify
  //! @return Classification result
  Standard_EXPORT Result SiDans(const gp_Pnt2d& thePoint) const;

  //! Classifies a point with explicit ON tolerance.
  //!
  //! Similar to SiDans() but uses the specified tolerance for boundary detection
  //! instead of the tolerances specified at construction.
  //!
  //! @param[in] thePoint The 2D point to classify
  //! @param[in] theTol   Tolerance for boundary detection
  //! @return Classification result
  Standard_EXPORT Result SiDans_OnMode(const gp_Pnt2d& thePoint, double theTol) const;

private:
  //! Internal classification in normalized coordinates.
  //!
  //! Performs point-in-polygon test using ray-casting algorithm.
  //! The point coordinates should be in the normalized [0,1] domain.
  //!
  //! @param[in] theX X coordinate in normalized space
  //! @param[in] theY Y coordinate in normalized space
  //! @return true if inside, false if outside
  bool internalSiDans(double theX, double theY) const;

  //! Internal classification with ON detection.
  //!
  //! Same as internalSiDans() but also detects if the point lies on the boundary.
  //!
  //! @param[in] theX X coordinate in normalized space
  //! @param[in] theY Y coordinate in normalized space
  //! @return Classification result
  Result internalSiDansOuOn(double theX, double theY) const;

  //! Builds the bucket index of the polygon's edges by their range in V (large polygons only).
  //!
  //! Whatever an edge can contribute to a classification at height Y (a crossing of the ray, a
  //! vertex or an edge within the tolerance) needs Y within the edge's V range widened by the
  //! tolerance. The edges are listed under every bucket their widened range meets, so a query
  //! reads one bucket instead of the whole polygon and finds the same edges, in the same order.
  void buildIndex();

  //! The edges that can matter at height theY: a list of edge numbers (1..myPointsCount, edge
  //! k runs from vertex k-1 to vertex k), or no list for every edge in turn.
  //! @param[out] theList  edge numbers, or nullptr for 1..theCount
  //! @param[out] theCount number of edges to visit
  void candidates(double theY, const int*& theList, int& theCount) const;

  //! Initializes the classifier with polygon data.
  //! @tparam TCol_Containers2d Container type (Array1 or Sequence)
  template <class TCol_Containers2d>
  void init(const TCol_Containers2d& thePnts2d,
            double                   theTolU,
            double                   theTolV,
            double                   theUMin,
            double                   theVMin,
            double                   theUMax,
            double                   theVMax);

  //! Copy constructor is deleted.
  CSLib_Class2d(const CSLib_Class2d&) = delete;

  //! Copy assignment operator is deleted.
  CSLib_Class2d& operator=(const CSLib_Class2d&) = delete;

private:
  NCollection_Array1<double> myPnts2dX;           //!< X coordinates (normalized)
  NCollection_Array1<double> myPnts2dY;           //!< Y coordinates (normalized)
  NCollection_Array1<int>    myBucketStart;       //!< First entry of each V bucket in myBucketEdges
  NCollection_Array1<int>    myBucketEdges;       //!< Edge numbers by V bucket, ascending in each
  int                        myBucketCount = 0;   //!< Number of V buckets, 0 for no index
  double                     myTolU        = 0.0; //!< Tolerance in U direction (normalized)
  double                     myTolV        = 0.0; //!< Tolerance in V direction (normalized)
  int                        myPointsCount = 0;   //!< Number of polygon vertices
  double                     myUMin        = 0.0; //!< Original minimum U bound
  double                     myVMin        = 0.0; //!< Original minimum V bound
  double                     myUMax        = 0.0; //!< Original maximum U bound
  double                     myVMax        = 0.0; //!< Original maximum V bound
};

#endif // _CSLib_Class2d_HeaderFile
