#
# CMakeLists.txt
#   setup

# List of files to compile 
set(SETUP_FILES)

# List of source files
set(SETUP_SOURCE_FILES)

# Process target name
set(SETUP_TARGET_NAME ProcessTruchasSetupFiles)

# 
set(SETUP_FPP_FLAGS 
          -I${TruchasExe_SOURCE_DIR}/utilities
	  -I${TruchasExe_SOURCE_DIR}/setup/restart
	  -I${TruchasExe_SOURCE_DIR}/setup/scalars)

# - base_types
set(SETUP_BASE_FILES
          setup/base_types/base_types_A_module.F90
          setup/base_types/matl_module.F90
          setup/base_types/matl_utilities.F90
          setup/base_types/parallel_scope.F90
          setup/base_types/var_vector_types.F90
          setup/base_types/zone_module.F90)
list(APPEND SETUP_FILES ${SETUP_BASE_FILES})	

# - bc 	
set(SETUP_BC_FILES
          setup/bc/bc_module.F90
          setup/bc/bc_data_module.F90
          setup/bc/bc_flag_module.F90
          setup/bc/bc_info_module.F90
          setup/bc/bc_kind_module.F90
          setup/bc/bc_set_module.F90
          setup/bc/bc_type_module.F90
          setup/bc/velocity_boundary_data_type.F90)
list(APPEND SETUP_FILES ${SETUP_BC_FILES})	

# - initialize
set(SETUP_INIT_FILES
          setup/initialize/init_module.F90)
list(APPEND SETUP_FILES ${SETUP_INIT_FILES})	

# - restart
set(SETUP_RESTART_FILES
          setup/restart/restart_driver.F90
          setup/restart/restart_utilities.F90
          setup/restart/restart_variables.F90)
list(APPEND SETUP_FILES ${SETUP_RESTART_FILES})	


# - scalars
set(SETUP_SCALARS_FILES
          setup/scalars/constants_module.F90
          setup/scalars/cutoffs_module.F90
          setup/scalars/debug_control_data.F90
          setup/scalars/parameter_module.F90)
list(APPEND SETUP_FILES ${SETUP_SCALARS_FILES})	

# - vof
set(SETUP_VOF_FILES
          setup/vof/interfaces_module.F90
          setup/vof/body_identifier_type.F90
          setup/vof/body_volume_initialize_routine.F90
          setup/vof/body_class.F90
          setup/vof/background_body_type.F90
          setup/vof/element_block_body_type.F90
          setup/vof/ellipsoid_body_type.F90
          setup/vof/plane_body_type.F90
          setup/vof/box_body_type.F90
          setup/vof/sphere_body_type.F90
          setup/vof/cylinder_body_type.F90
          setup/vof/body_factories.F90
          setup/vof/body_namelist.F90
          )
list(APPEND SETUP_FILES ${SETUP_VOF_FILES})	


# - miscellenous files
list(APPEND SETUP_FILES
            setup/setup_module.F90
            setup/random_module.F90
	    setup/overwrite/overwrite_module.F90)

# Process Fortran files
fortran_preprocess_files(SETUP_SOURCE_FILES
                         FILES ${SETUP_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${SETUP_FPP_FLAGS}
			 PROCESS_TARGET ${SETUP_TARGET_NAME})

include(BuildWhitespaceString)
set(fc_flags -I${PGSLib_MODULE_DIR})
list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
build_whitespace_string(SETUP_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${SETUP_SOURCE_FILES} PROPERTIES
                              COMPILE_FLAGS ${SETUP_COMPILE_FLAGS})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SETUP_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS ${SETUP_TARGET_NAME})
