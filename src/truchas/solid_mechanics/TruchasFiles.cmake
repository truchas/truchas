# Truchas files in directory
#   solid_mechanics

set(SM_SOURCE_FILES
  solid_mechanics/solid_mechanics_type.F90
  solid_mechanics/integration_geometry_type.F90
  solid_mechanics/integration_cell_type.F90
  )

# Define compile flags
include(BuildWhitespaceString)
set(sm_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND sm_source_flags "-assume realloc_lhs")
endif()

build_whitespace_string(sm_source_flags_str ${sm_source_flags})
set_source_files_properties(${SM_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${sm_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SM_SOURCE_FILES})
