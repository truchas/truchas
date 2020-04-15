# Truchas files in directory
#   functions

set(FUNC_SOURCE_FILES
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
  list(APPEND FUNC_SOURCE_FILES functions/dl_scalar_func_type.F90)
endif()


# Define compile flags
include(BuildWhitespaceString)
set(func_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND func_source_flags "-standard-semantics -assume nostd_mod_proc_name")
endif()

build_whitespace_string(func_source_flags_str ${func_source_flags})
set_source_files_properties(${FUNC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${func_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FUNC_SOURCE_FILES})
