# Truchas files in directory
#   functions

set(FLOW_SOURCE_FILES
        flow/flow_domain_types.F90
        flow/flow_input_utils.F90
        flow/fischer_guess_type.F90
	flow/flow_bc_factory_type.F90
	flow/flow_bc_type.F90
	flow/flow_driver.F90
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
        flow/flow_surface_tension_bc_type.F90
        flow/flow_solver_namelists.F90
        flow/turbulence_namelist.F90
	)


# Define compile flags
include(BuildWhitespaceString)
set(flow_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})

build_whitespace_string(flow_source_flags_str ${flow_source_flags})
set_source_files_properties(${FLOW_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${flow_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${FLOW_SOURCE_FILES})
