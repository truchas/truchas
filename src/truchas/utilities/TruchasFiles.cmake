# Truchas files in directory
#   utilities

# List of files to  process
set(UTIL_FILES)

# List of files to add to the Truchas library
set(UTIL_SOURCE_FILES)


set(UTIL_FILES
         utilities/ArrayAllocate_Module.F90
         utilities/f90_assert.F90
         utilities/file_utility.F90
         utilities/graph_type.F90
         utilities/input_utilities.F90
         utilities/integer_set_type.F90
         utilities/kinds.F90
         utilities/linear_module.F90
         utilities/lnorm_module.F90
         utilities/lu_solve_module.F90
         utilities/old_mesh_gmv.F90
         utilities/permutations.F90
         utilities/process_info_module.F90
         utilities/sort_module.F90
         utilities/sort_utilities.F90
         utilities/string_utilities.F90
         utilities/tabular_utilities.F90
         utilities/tensor_module.F90
         utilities/timing_tree.F90
         utilities/truchas_env.F90
         utilities/truchas_logging_services.F90
         utilities/truchas_timing.F90
         utilities/utilities_module.F90
         utilities/var_vector_module.F90
	 utilities/gmv/fgmvwrite.F90)

set(UTIL_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(UTIL_SOURCE_FILES
                         FILES ${UTIL_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${UTIL_FPP_FLAGS}
			 PROCESS_TARGET ProcessTruchasUtilFiles)
set_source_files_properties(${UTIL_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${PGSLib_MODULE_DIR})

# Add the C files 
list(APPEND UTIL_SOURCE_FILES
            utilities/get_process_size.c
            utilities/make_directory.c
            utilities/gmv/fgmvwrite_c.c
            utilities/gmv/gmvwrite.c)
set(gmv_compile_flags "-I${Truchas_FCIface_INCLUDE_DIR} -I${TruchasExe_SOURCE_DIR}/utilities/gmv")	  
set_source_files_properties(utilities/gmv/fgmvwrite_c.c PROPERTIES
                            COMPILE_FLAGS ${gmv_compile_flags})
set_source_files_properties(utilities/make_directory.c PROPERTIES
                            COMPILE_FLAGS -I${Truchas_FCIface_INCLUDE_DIR})

include(CheckFunctionExists)
check_function_exists(getpid HAVE_GETPID)
if (HAVE_GETPID)
  set_source_files_properties(utilities/get_process_size.c PROPERTIES
                              COMPILE_FLAGS -DLINUX)
endif(HAVE_GETPID)			    

# Update the Truchas library file list and targets		       
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${UTIL_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ProcessTruchasUtilFiles)

