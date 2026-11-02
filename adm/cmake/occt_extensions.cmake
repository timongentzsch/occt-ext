# Options of the extension packages: TKAlgoExt, and the packages GeomAlgoExt (in TKGeomAlgo) and
# BRepAlgoExt (in TKMesh). Included after the toolkits have been added.

option(BUILD_EXT_KEEP_EXCEPTIONS
       "Compile the extension packages with C++ exceptions and assertions on, whatever OCCT's own flags say (-DNo_Exception, -DNDEBUG)"
       ON)
set(BREPBLENDEXT_POLICY_HEADER "" CACHE FILEPATH
    "Header that defines the section laws of the rolling-ball blend (see BRepBlendExt_Policy.hxx); force-included into the BRepBlendExt sources. Empty: the built-in policy")
set(BREPBLENDEXT_POLICY_CLASS "BRepExt::blend2::NativePolicy" CACHE STRING
    "The class of that header that BRepBlendExt uses")

if (BUILD_EXT_KEEP_EXCEPTIONS)
  if (MSVC)
    set (EXT_FLAGS /UNo_Exception /UNDEBUG)
  else()
    set (EXT_FLAGS -UNo_Exception -UNDEBUG)
  endif()
  if (TARGET TKAlgoExt)
    target_compile_options(TKAlgoExt PRIVATE ${EXT_FLAGS})
  endif()
  # Compiled into toolkits that otherwise keep OCCT's flags.
  file(GLOB EXT_PACKAGE_SOURCES
       "${OCCT_ROOT_DIR}/src/ModelingAlgorithms/TKGeomAlgo/GeomAlgoExt/*.cxx"
       "${OCCT_ROOT_DIR}/src/ModelingAlgorithms/TKMesh/BRepAlgoExt/*.cxx")
  foreach (EXT_TOOLKIT TKGeomAlgo TKMesh)
    if (TARGET ${EXT_TOOLKIT})
      set_source_files_properties(${EXT_PACKAGE_SOURCES} TARGET_DIRECTORY ${EXT_TOOLKIT}
                                  PROPERTIES COMPILE_OPTIONS "${EXT_FLAGS}")
    endif()
  endforeach()
endif()

if (BREPBLENDEXT_POLICY_HEADER AND TARGET TKAlgoExt)
  target_compile_definitions(TKAlgoExt PRIVATE BREPBLENDEXT_POLICY=${BREPBLENDEXT_POLICY_CLASS})
  file(GLOB BLEND_SOURCES "${OCCT_ROOT_DIR}/src/ModelingAlgorithms/TKAlgoExt/BRepBlendExt/*.cxx")
  if (MSVC)
    set (BLEND_INCLUDE "/FI${BREPBLENDEXT_POLICY_HEADER}")
  else()
    set (BLEND_INCLUDE "-include;${BREPBLENDEXT_POLICY_HEADER}")
  endif()
  set_source_files_properties(${BLEND_SOURCES} TARGET_DIRECTORY TKAlgoExt
                              PROPERTIES COMPILE_OPTIONS "${BLEND_INCLUDE}")
endif()
