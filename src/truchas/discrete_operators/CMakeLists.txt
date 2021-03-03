target_sources(truchas PRIVATE
  discrete_derivatives.F90
  discrete_op_module.F90
  discrete_ops_data.F90
  do_base_types.F90
  do_discrete_operators.F90
  do_interface.F90
  do_solve_module.F90
  do_solve_specifier.F90
  do_update_module.F90
  ff_discrete_ops_data.F90
  cell_grad_type.F90
)

# Define compile flags
#include(BuildWhitespaceString)
#set(DISOP_COMPILE_FLAGS
#  -I${PGSLib_MODULE_DIR}
#  -I${PETACA_MODULE_DIR}
#  -I${Truchas_utilities_dir})

#build_whitespace_string(DISOP_COMPILE_FLAGS_STR ${DISOP_COMPILE_FLAGS})
#set_source_files_properties(${DISOP_SOURCE_FILES} PROPERTIES
#                            COMPILE_FLAGS ${DISOP_COMPILE_FLAGS_STR})


#list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISOP_SOURCE_FILES})
