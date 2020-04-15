# Truchas files in directory
#   legacy_mesh

set(LEGACYMESH_SOURCE_FILES
    legacy_mesh/legacy_mesh_api.F90
    legacy_mesh/common_impl.F90
    legacy_mesh/mesh_impl.F90
    legacy_mesh/vertex_impl.F90
    legacy_mesh/cell_impl.F90
    legacy_mesh/ee_gather_impl.F90
    legacy_mesh/en_gather_impl.F90
    legacy_mesh/nn_gather_impl.F90
    legacy_mesh/mesh_face_set_impl.F90
    legacy_mesh/legacy_geometry.F90
    )

# Define compile flags
set(LEGACYMESH_COMPILE_FLAGS
  -I${PGSLib_MODULE_DIR}
  -I${Truchas_utilities_dir})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    list(APPEND LEGACYMESH_COMPILE_FLAGS "-assume realloc_lhs")
endif()
include(BuildWhitespaceString)
build_whitespace_string(LEGACYMESH_COMPILE_FLAGS_STR ${LEGACYMESH_COMPILE_FLAGS})
set_source_files_properties(${LEGACYMESH_SOURCE_FILES} PROPERTIES
    COMPILE_FLAGS ${LEGACYMESH_COMPILE_FLAGS_STR})
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${LEGACYMESH_SOURCE_FILES})
