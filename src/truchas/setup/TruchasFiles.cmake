#
# CMakeLists.txt
#   setup

# - base_types
set(SETUP_SOURCE_FILES
          setup/base_types/base_types_A_module.F90
          setup/base_types/matl_module.F90
          setup/base_types/matl_utilities.F90
          setup/base_types/parallel_scope.F90
          setup/base_types/var_vector_types.F90
          setup/base_types/zone_module.F90)

# - bc
list(APPEND SETUP_SOURCE_FILES
          setup/bc/bc_module.F90
          setup/bc/bc_data_module.F90
          setup/bc/bc_flag_module.F90
          setup/bc/bc_info_module.F90
          setup/bc/bc_kind_module.F90
          setup/bc/bc_set_module.F90
          setup/bc/bc_type_module.F90
          setup/bc/velocity_boundary_data_type.F90)

# - initialize
list(APPEND SETUP_SOURCE_FILES
          setup/initialize/init_module.F90)

# - restart
list(APPEND SETUP_SOURCE_FILES
          setup/restart/restart_driver.F90
          setup/restart/restart_utilities.F90
          setup/restart/restart_variables.F90)

# - scalars
list(APPEND SETUP_SOURCE_FILES
          setup/scalars/constants_module.F90
          setup/scalars/cutoffs_module.F90
          setup/scalars/debug_control_data.F90
          setup/scalars/parameter_module.F90)

# - vof
list(APPEND SETUP_SOURCE_FILES
          setup/vof/interfaces_module.F90
          setup/vof/body_identifier_type.F90
          setup/vof/compute_body_volumes_proc.F90
          setup/vof/body_class.F90
          setup/vof/background_body_type.F90
          setup/vof/element_block_body_type.F90
          setup/vof/ellipsoid_body_type.F90
          setup/vof/ellipse_body_type.F90
          setup/vof/plane_body_type.F90
          setup/vof/box_body_type.F90
          setup/vof/sphere_body_type.F90
          setup/vof/cylinder_body_type.F90
          setup/vof/body_factories.F90
          setup/vof/body_namelist.F90
          )


# - miscellenous files
list(APPEND SETUP_SOURCE_FILES
            setup/setup_module.F90
            setup/random_module.F90
	    setup/overwrite/overwrite_module.F90)


include(BuildWhitespaceString)
set(fc_flags
  -I${PGSLib_MODULE_DIR}
  -I${Truchas_utilities_dir})
build_whitespace_string(SETUP_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${SETUP_SOURCE_FILES} PROPERTIES
                              COMPILE_FLAGS ${SETUP_COMPILE_FLAGS})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SETUP_SOURCE_FILES})
