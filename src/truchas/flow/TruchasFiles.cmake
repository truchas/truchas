# Truchas files in directory
#   functions

# List of files to  process
set(FLOW_FILES)

# List of files to add to the Truchas library
set(FLOW_SOURCE_FILES)

# Process target name
set(FLOW_TARGET_NAME ProcessTruchasFlowFiles)

set(FLOW_FILES
        flow/flow_domain_types.F90
        flow/flow_input_utils.F90
        flow/fischer_guess_type.F90
	flow/flow_bc_factory_type.F90
	flow/flow_bc_type.F90
	flow/flow_driver.F90
	flow/flow_mesh_type.F90
	flow/flow_operators.F90
	flow/turbulence_model_class.F90
	flow/default_turb_model_type.F90
	flow/algebraic_turb_model_type.F90
	flow/turbulence_models.F90
	flow/flow_projection_type.F90
	flow/flow_prediction_type.F90
	flow/flow_props_type.F90
	flow/flow_type.F90
        flow/flow_namelist.F90
        flow/flow_bc_namelist.F90
        flow/flow_predictor_namelist.F90
        flow/flow_corrector_namelist.F90
        flow/turbulence_namelist.F90
	)

set(FLOW_FPP_FLAGS
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(FLOW_SOURCE_FILES
                         FILES ${FLOW_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${FLOW_FPP_FLAGS}
                         PROCESS_TARGET ${FLOW_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(flow_source_flags -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND flow_source_flags "-assume realloc_lhs")
endif()

build_whitespace_string(flow_source_flags_str ${flow_source_flags})
set_source_files_properties(${FLOW_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${flow_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FLOW_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${FLOW_TARGET_NAME})
