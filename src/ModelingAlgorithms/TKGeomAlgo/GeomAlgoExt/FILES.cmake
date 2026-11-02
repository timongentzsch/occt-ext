set(OCCT_GeomAlgoExt_FILES_LOCATION "${CMAKE_CURRENT_LIST_DIR}")
set(OCCT_GeomAlgoExt_FILES
  GeomAlgoExt_PatchFit.hxx
  GeomAlgoExt_PatchFit.cxx
  GeomAlgoExt_PatchRim.hxx
  GeomAlgoExt_PatchRim.cxx
  GeomAlgoExt_PatchCompact.hxx
  GeomAlgoExt_PatchCompact.cxx
  GeomAlgoExt_Failure.hxx
  GeomAlgoExt_Failure.cxx
  GeomAlgoExt_Cancel.hxx
  GeomAlgoExt_Phase.hxx
  GeomAlgoExt_FitDiagnostics.hxx
  GeomAlgoExt_FitDiagnostics.cxx
  GeomAlgoExt_PlanarCurve.hxx
  GeomAlgoExt_SurfaceJet.hxx
  GeomAlgoExt_Resolution.hxx
  GeomAlgoExt_Coordinates.hxx
  GeomAlgoExt_Progress.hxx
  GeomAlgoExt_Band.hxx
  GeomAlgoExt_Pair.hxx
  GeomAlgoExt_Bernstein.hxx
  GeomAlgoExt_CubicThrough.hxx
  GeomAlgoExt_Hermite.hxx
  GeomAlgoExt_Interval.hxx
  GeomAlgoExt_ParameterRectangle.hxx
  GeomAlgoExt_SurfaceBounds.hxx
  GeomAlgoExt_UniformKnots.hxx
  GeomAlgoExt_ModelTolerance.hxx
  GeomAlgoExt_Timing.hxx
)

# Native defaults must not resolve a WebAssembly host's imported hook symbols.
# WebAssembly applications supply their own direct cancellation/observer hooks.
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  list(APPEND OCCT_GeomAlgoExt_FILES GeomAlgoExt_Progress.cxx)
endif()
