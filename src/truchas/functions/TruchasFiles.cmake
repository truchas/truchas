# Truchas files in directory
#   functions

# List of files to  process
set(FUNC_FILES)

# List of files to add to the Truchas library
set(FUNC_SOURCE_FILES)

# Process target name
set(FUNC_TARGET_NAME ProcessTruchasFunctionFiles)


set(FUNC_FILES
        functions/scalar_func_table.F90
	functions/vector_func_table.F90
        functions/function_namelist.F90
	functions/vfunction_namelist.F90
        functions/scalar_func_class.F90
        functions/const_scalar_func_type.F90
        functions/poly_scalar_func_type.F90
        functions/mpoly_scalar_func_type.F90
        functions/tabular_scalar_func_type.F90
        functions/fptr_scalar_func_type.F90
        functions/scalar_func_factories.F90
        functions/scalar_func_containers.F90
	functions/vector_func_containers.F90
        functions/scalar_func_tools.F90
        functions/scalar_func_map_type.F90
	functions/vector_func_map_type.F90
        functions/vector_func_class.F90
        functions/const_vector_func_type.F90
        functions/tabular_vector_func_type.F90
        functions/fptr_vector_func_type.F90
        functions/vector_func_factories.F90
        )

if(ENABLE_DYNAMIC_LOADING)
  list(APPEND FUNC_FILES functions/dl_scalar_func_type.F90)
endif()

set(FUNC_FPP_FLAGS
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(FUNC_SOURCE_FILES
                         FILES ${FUNC_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${FUNC_FPP_FLAGS}
                         PROCESS_TARGET ${FUNC_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(func_source_flags -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND func_source_flags "-standard-semantics -assume nostd_mod_proc_name")
endif()

build_whitespace_string(func_source_flags_str ${func_source_flags})
set_source_files_properties(${FUNC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${func_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FUNC_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${FUNC_TARGET_NAME})
