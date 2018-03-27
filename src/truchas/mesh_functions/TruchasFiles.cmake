# Truchas files in directory
#   mesh_functions

# List of files to  process
set(MESH_FUNC_FILES)

# List of files to add to the Truchas library
set(MESH_FUNC_SOURCE_FILES)

# Process target name
set(MESH_FUNC_TARGET_NAME ProcessTruchasMeshFunctionFiles)

set(MESH_FUNC_FILES
    mesh_functions/bndry_func_class.F90
    mesh_functions/bndry_face_func_type.F90
    mesh_functions/bndry_face_group_builder_type.F90
)

set(MESH_FUNC_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(MESH_FUNC_SOURCE_FILES
                         FILES ${MESH_FUNC_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${MESH_FUNC_FPP_FLAGS}
                         PROCESS_TARGET ${MESH_FUNC_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(mesh_func_source_flags -I${PGSLib_MODULE_DIR} -I${Danu_Fortran_MODULE_DIR} -I${PETACA_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND mesh_func_source_flags "-standard-semantics -assume nostd_mod_proc_name")
endif()

build_whitespace_string(mesh_func_source_flags_str ${mesh_func_source_flags})
set_source_files_properties(${MESH_FUNC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${mesh_func_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${MESH_FUNC_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${MESH_FUNC_TARGET_NAME})


