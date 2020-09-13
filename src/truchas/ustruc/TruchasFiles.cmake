# Truchas files in directory
#   functions

set(USTRUC_SOURCE_FILES
        ustruc/ustruc_driver.F90
        ustruc/ustruc_model_type.F90
        ustruc/ustruc_comp_class.F90
        ustruc/ustruc_core_type.F90
        ustruc/ustruc_plugin_class.F90
        ustruc/ustruc_vel1_type.F90
        ustruc/ustruc_time_type.F90
        ustruc/ustruc_gv0_type.F90
        ustruc/ustruc_gv1_type.F90
        ustruc/ustruc_comp_factory.F90
        ustruc/serialization_tools.F90
        )


# Define compile flags
include(BuildWhitespaceString)
set(ustruc_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})

build_whitespace_string(ustruc_source_flags_str ${ustruc_source_flags})
set_source_files_properties(${USTRUC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${ustruc_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${USTRUC_SOURCE_FILES})
