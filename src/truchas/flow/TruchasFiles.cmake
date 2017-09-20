# Truchas files in directory
#   functions

# List of files to  process
set(FLOW_FILES)

# List of files to add to the Truchas library
set(FLOW_SOURCE_FILES)

# Process target name
set(FLOW_TARGET_NAME ProcessTruchasFlowFiles)

set(FLOW_FILES
        flow/flow_driver.F90
        flow/flow_model.F90
        flow/vof_model.F90
        # flow/flow_core_type.F90
        # flow/flow_plugin_class.F90
        # flow/flow_vel1_type.F90
        # flow/flow_time_type.F90
        # flow/flow_gv0_type.F90
        # flow/flow_gv1_type.F90
        # flow/flow_comp_factory.F90
        # flow/serialization_tools.F90
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
