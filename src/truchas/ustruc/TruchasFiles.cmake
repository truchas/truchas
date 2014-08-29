# Truchas files in directory
#   functions

# List of files to  process
set(USTRUC_FILES)

# List of files to add to the Truchas library
set(USTRUC_SOURCE_FILES)

# Process target name
set(USTRUC_TARGET_NAME ProcessTruchasUstrucFiles)

set(USTRUC_FILES
        ustruc/ustruc_driver.F90
        ustruc/ustruc_model_type.F90
        ustruc/ustruc_comp_class.F90
        ustruc/ustruc_core_type.F90
        ustruc/ustruc_plugin_class.F90
        ustruc/ustruc_vel1_type.F90
        ustruc/ustruc_time_type.F90
        ustruc/ustruc_gv1_type.F90
        ustruc/ustruc_comp_factory.F90
        )

set(USTRUC_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(USTRUC_SOURCE_FILES
                         FILES ${USTRUC_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${USTRUC_FPP_FLAGS}
                         PROCESS_TARGET ${USTRUC_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(ustruc_source_flags -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR})
if(Fortran_COMPILER_IS_INTEL)
  list(APPEND ustruc_source_flags "-assume realloc_lhs")
endif()

build_whitespace_string(ustruc_source_flags_str ${ustruc_source_flags})
set_source_files_properties(${USTRUC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${ustruc_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${USTRUC_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${USTRUC_TARGET_NAME})


