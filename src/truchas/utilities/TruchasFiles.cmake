# Truchas files in directory
#   utilities

set(UTIL_SOURCE_FILES
         utilities/f90_assert.F90
         utilities/file_utility.F90
         utilities/graph_type.F90
         utilities/input_utilities.F90
         utilities/integer_set_type.F90
         utilities/integer_set_type_wavl.F90
         utilities/string_set_type.F90
         utilities/kinds.F90
         utilities/lnorm_module.F90
         utilities/lu_solve_module.F90
         utilities/old_mesh_gmv.F90
         utilities/permutations.F90
         utilities/process_info_module.F90
         utilities/signal_handler.F90
         utilities/sort_module.F90
         utilities/sort_utilities.F90
         utilities/string_utilities.F90
         utilities/tabular_utilities.F90
         utilities/tensor_module.F90
         utilities/truchas_env.F90
         utilities/truchas_logging_services.F90
         utilities/truchas_timers.F90
         utilities/utilities_module.F90
         utilities/var_vector_module.F90
         utilities/gmv/gmvwrite_c_binding.F90
	 utilities/f08_intrinsics.F90)

set_source_files_properties(${UTIL_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLib_MODULE_DIR}")

# Add the C files
list(APPEND UTIL_SOURCE_FILES
            utilities/get_process_size.c
            utilities/make_directory.c
            utilities/gmv/gmvwrite.c)

if(CMAKE_SYSTEM_NAME MATCHES Linux)
  set_source_files_properties(utilities/get_process_size.c PROPERTIES
                              COMPILE_DEFINITIONS LINUX)
endif()

# Update the Truchas library file list and targets
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${UTIL_SOURCE_FILES})
