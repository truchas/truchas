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
          setup/base_types/base_types_module.F90
          setup/base_types/matl_module.F90
          setup/base_types/matl_utilities.F90
          setup/base_types/mesh_module.F90
          setup/base_types/parallel_scope.F90
          setup/base_types/probe_module.F90
          setup/base_types/tempGrad_module.F90
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
          setup/bc/bc_type_module.F90)
list(APPEND SETUP_FILES ${SETUP_BC_FILES})	

# - initialize
set(SETUP_INIT_FILES
          setup/initialize/init_module.F90)
list(APPEND SETUP_FILES ${SETUP_INIT_FILES})	

# - mesh
set(SETUP_MESH_FILES
          setup/mesh/mesh_decomposition_module.F90
          setup/mesh/mesh_distribute_module.F90
          setup/mesh/mesh_gen_data.F90
          setup/mesh/mesh_gen_module.F90
          setup/mesh/mesh_partition_module.F90
          setup/mesh/mesh_quality_module.F90
          setup/mesh/mesh_tests.F90
          setup/mesh/mesh_utilities.F90
          setup/mesh/partitioner_data.F90
          setup/mesh/two_level_partition.F90)
list(APPEND SETUP_FILES ${SETUP_MESH_FILES})	

# - restart
set(SETUP_RESTART_FILES
          setup/restart/restart_driver.F90
          setup/restart/restart_utilities.F90
          setup/restart/restart_variables.F90)
list(APPEND SETUP_FILES ${SETUP_RESTART_FILES})	


# - scalars
set(SETUP_SCALARS_FILES
          setup/scalars/code_module.F90
          setup/scalars/constants_module.F90
          setup/scalars/cutoffs_module.F90
          setup/scalars/debug_control_data.F90
          setup/scalars/parameter_module.F90
          setup/scalars/scalars_module.F90)
list(APPEND SETUP_FILES ${SETUP_SCALARS_FILES})	

# - vof
set(SETUP_VOF_FILES
          setup/vof/interfaces_module.F90
          setup/vof/tally_module.F90
          setup/vof/vof_init.F90)         
list(APPEND SETUP_FILES ${SETUP_VOF_FILES})	


# - miscellenous files
list(APPEND SETUP_FILES
            setup/setup_module.F90
            setup/random_module.F90
	    setup/cell_geometry/cell_geometry_module.F90
	    setup/overwrite/overwrite_module.F90)

# Process Fortran files
fortran_preprocess_files(SETUP_SOURCE_FILES
                         FILES ${SETUP_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${SETUP_FPP_FLAGS}
			 PROCESS_TARGET ${SETUP_TARGET_NAME})

include(BuildWhitespaceString)
if(ENABLE_PGSLib)
  set(fc_flags -I${PGSLib_MODULE_DIR})
endif()
if(ENABLE_Danu)
  list(APPEND fc_flags -I${Danu_Fortran_MODULE_DIR})
endif()
if(ENABLE_UbikSolve)
  list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
endif()
list(APPEND fc_flags -I${NETCDF_INCLUDE_DIR})
build_whitespace_string(SETUP_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${SETUP_SOURCE_FILES} PROPERTIES
                              COMPILE_FLAGS ${SETUP_COMPILE_FLAGS})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SETUP_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS ${SETUP_TARGET_NAME})

# Add C files to Truchas source file list
list(APPEND Truchas_LIBRARY_SOURCE_FILES
                    setup/mesh/chaco_f90_wrapper.c)
set(chaco_cflags "-I${Truchas_FCIface_INCLUDE_DIR}")
if(ENABLE_Chaco)
  set(chaco_cflags "${chaco_cflags} -DUSE_CHACO")
endif()  
set_source_files_properties(setup/mesh/chaco_f90_wrapper.c PROPERTIES
                            COMPILE_FLAGS ${chaco_cflags})



