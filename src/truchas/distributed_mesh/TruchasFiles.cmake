# Truchas files in directory
#   distributed_mesh

set(DISMESH_SOURCE_FILES
           distributed_mesh/bitfield_type.F90
           distributed_mesh/base_mesh_class.F90
           distributed_mesh/unstr_base_mesh_class.F90
           distributed_mesh/facet_hash_type.F90
           distributed_mesh/facet_table_type.F90
           distributed_mesh/face_neighbor_table_type.F90
           distributed_mesh/index_partitioning.F90
           distributed_mesh/mesh_manager.F90
           distributed_mesh/parallel_communication.F90
           distributed_mesh/parallel_permutations.F90
           distributed_mesh/cell_topology.F90
           distributed_mesh/cell_geometry.F90
           distributed_mesh/unstr_mesh_type.F90
           distributed_mesh/unstr_mesh_factory.F90
           distributed_mesh/unstr_mesh_tools.F90
           distributed_mesh/unstr_mesh_gmv.F90
           distributed_mesh/ext_exodus_mesh_type.F90
           distributed_mesh/exodus_mesh_tools.F90
           distributed_mesh/simplex_topology.F90
           distributed_mesh/simplex_geometry.F90
	   distributed_mesh/simpl_mesh_type.F90
           distributed_mesh/simpl_mesh_factory.F90
           distributed_mesh/simpl_mesh_tools.F90
           distributed_mesh/simpl_mesh_gmv.F90
           )

# Define compile flags
set(DISMESH_COMPILE_FLAGS -I${PGSLib_MODULE_DIR}
        -I${TruchasExe_SOURCE_DIR}/distributed_mesh
	-I${TruchasExe_SOURCE_DIR}/utilities
	)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND DISMESH_COMPILE_FLAGS "-assume realloc_lhs")
endif()
include(BuildWhitespaceString)
build_whitespace_string(DISMESH_COMPILE_FLAGS_STR ${DISMESH_COMPILE_FLAGS})
set_source_files_properties(${DISMESH_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${DISMESH_COMPILE_FLAGS_STR})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISMESH_SOURCE_FILES})
