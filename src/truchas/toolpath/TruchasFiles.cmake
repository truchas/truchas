# Truchas files in directory
#   toolpath

# List of files to  process
set(TOOLPATH_FILES)

# List of files to add to the Truchas library
set(TOOLPATH_SOURCE_FILES)

# Process target name
set(TOOLPATH_TARGET_NAME ProcessTruchasToolpathFiles)

set(TOOLPATH_FILES
        toolpath/toolpath_type.F90
        toolpath/xyz_motion_class.F90
        toolpath/dwell_xyz_motion_type.F90
        toolpath/linear_xyz_motion_type.F90
        toolpath/toolpath_factory.F90
        toolpath/toolpath_namelist.F90
        toolpath/toolpath_table.F90
        )

set(TOOLPATH_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(TOOLPATH_SOURCE_FILES
                         FILES ${TOOLPATH_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${TOOLPATH_FPP_FLAGS}
                         PROCESS_TARGET ${TOOLPATH_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(toolpath_source_flags -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND toolpath_source_flags "-assume realloc_lhs")
endif()

build_whitespace_string(toolpath_source_flags_str ${toolpath_source_flags})
set_source_files_properties(${TOOLPATH_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${toolpath_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${TOOLPATH_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${TOOLPATH_TARGET_NAME})


