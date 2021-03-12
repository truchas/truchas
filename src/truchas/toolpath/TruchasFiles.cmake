# Truchas files in directory
#   toolpath

set(TOOLPATH_SOURCE_FILES
        toolpath/toolpath_type.F90
        toolpath/xyz_motion_class.F90
        toolpath/dwell_xyz_motion_type.F90
        toolpath/linear_xyz_motion_type.F90
        toolpath/toolpath_factory.F90
        toolpath/toolpath_namelist.F90
        toolpath/toolpath_table.F90
        )


# Define compile flags
include(BuildWhitespaceString)
set(toolpath_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})

build_whitespace_string(toolpath_source_flags_str ${toolpath_source_flags})
set_source_files_properties(${TOOLPATH_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${toolpath_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${TOOLPATH_SOURCE_FILES})
