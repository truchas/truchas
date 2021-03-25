# Truchas files in directory
#   solid_mechanics

set(SM_SOURCE_FILES
  solid_mechanics/solid_mechanics_driver.F90
  solid_mechanics/solid_mechanics_type.F90
  solid_mechanics/solid_mechanics_namelist.F90
  solid_mechanics/sm_bc_namelist.F90
  solid_mechanics/integration_geometry_type.F90
  solid_mechanics/integration_cell_type.F90
  solid_mechanics/bndry_ip_func_type.F90

  solid_mechanics/sm_bc_utilities.F90
  solid_mechanics/sm_normal_displacement_bc_type.F90
  solid_mechanics/sm_normal_traction_bc_type.F90
  solid_mechanics/sm_gap_contact_bc_type.F90
  solid_mechanics/sm_bc_type.F90
  solid_mechanics/sm_bc_node_types.F90
  solid_mechanics/sm_bc_face_type.F90
  #solid_mechanics/sm_bc_list_type.F90

  solid_mechanics/sm_model_type.F90
  solid_mechanics/sm_ds_precon_type.F90
  solid_mechanics/sm_nlsol_model_type.F90
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
