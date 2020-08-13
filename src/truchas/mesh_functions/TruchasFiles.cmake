# Truchas files in directory
#   mesh_functions

set(MESH_FUNC_SOURCE_FILES
    mesh_functions/bndry_vfunc_class.F90
    mesh_functions/bndry_face_vfunc_type.F90
    mesh_functions/bndry_func1_class.F90
    mesh_functions/bndry_func2_class.F90
    mesh_functions/intfc_func2_class.F90
    mesh_functions/bndry_face_func_type.F90
    #mesh_functions/bc_factory_type.F90
    mesh_functions/bndry_face_group_builder_type.F90
    mesh_functions/intfc_link_group_builder_type.F90
    mesh_functions/cell_group_builder_type.F90
    mesh_functions/htc_bndry_func_type.F90
    mesh_functions/rad_bndry_func_type.F90
    mesh_functions/htc_intfc_func_type.F90
    mesh_functions/rad_intfc_func_type.F90
    mesh_functions/scalar_mesh_func_class.F90
    mesh_functions/scalar_cell_func1_type.F90
    mesh_functions/scalar_cell_func2_type.F90
)


# Define compile flags
include(BuildWhitespaceString)
set(mesh_func_source_flags
  -I${PGSLib_MODULE_DIR} -I${PETACA_MODULE_DIR} -I${Truchas_utilities_dir})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND mesh_func_source_flags "-standard-semantics -assume nostd_mod_proc_name")
endif()

build_whitespace_string(mesh_func_source_flags_str ${mesh_func_source_flags})
set_source_files_properties(${MESH_FUNC_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${mesh_func_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${MESH_FUNC_SOURCE_FILES})
