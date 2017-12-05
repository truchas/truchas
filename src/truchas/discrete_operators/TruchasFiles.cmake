# Truchas files in directory
#   discrete_operators

# List of files to  process
set(DISOP_FILES)

# List of files to add to the Truchas library
set(DISOP_SOURCE_FILES)

# Process target name
set(DISOP_TARGET_NAME ProcessTruchasDisOpFiles)


set(DISOP_FILES
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
	 discrete_operators/flow_operators.F90)

set(DISOP_FPP_FLAGS
	-I${TruchasExe_SOURCE_DIR}/utilities ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(DISOP_SOURCE_FILES
                         FILES ${DISOP_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${DISOP_FPP_FLAGS}
			 PROCESS_TARGET ${DISOP_TARGET_NAME})
# Define compile flags
include(BuildWhitespaceString)
set(DISOP_COMPILE_FLAGS -I${PGSLib_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR}
                        -I${UbikSolve_MODULE_DIR} -I${PETACA_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND DISOP_COMPILE_FLAGS "-assume realloc_lhs") # for cell_grad_type.F90
endif()

build_whitespace_string(DISOP_COMPILE_FLAGS_STR ${DISOP_COMPILE_FLAGS})
set_source_files_properties(${DISOP_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${DISOP_COMPILE_FLAGS_STR})


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISOP_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS ${DISOP_TARGET_NAME})
