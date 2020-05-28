# Truchas files in directory
#   discrete_operators

set(DISOP_SOURCE_FILES
         discrete_operators/discrete_derivatives.F90
         discrete_operators/discrete_op_module.F90
         discrete_operators/discrete_ops_data.F90
         discrete_operators/do_base_types.F90
         discrete_operators/do_discrete_operators.F90
         discrete_operators/do_interface.F90
         discrete_operators/do_solve_module.F90
         discrete_operators/do_solve_specifier.F90
         discrete_operators/do_update_module.F90
         discrete_operators/ff_discrete_ops_data.F90
         discrete_operators/cell_grad_type.F90
	 )

# Define compile flags
include(BuildWhitespaceString)
set(DISOP_COMPILE_FLAGS
  -I${PGSLib_MODULE_DIR}
  -I${UbikSolve_MODULE_DIR}
  -I${PETACA_MODULE_DIR}
  -I${Truchas_utilities_dir})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND DISOP_COMPILE_FLAGS "-assume realloc_lhs") # for cell_grad_type.F90
endif()

build_whitespace_string(DISOP_COMPILE_FLAGS_STR ${DISOP_COMPILE_FLAGS})
set_source_files_properties(${DISOP_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${DISOP_COMPILE_FLAGS_STR})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISOP_SOURCE_FILES})
