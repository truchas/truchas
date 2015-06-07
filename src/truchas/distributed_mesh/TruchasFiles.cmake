# Truchas files in directory
#   distributed_mesh

# List of files to  process
set(DISMESH_FILES)

# List of files to add to the Truchas library
set(DISMESH_SOURCE_FILES)

# Process target name
set(DISMESH_TARGET_NAME ProcessTruchasDisMeshFiles)

set(DISMESH_FILES
           distributed_mesh/bitfield_type.F90
           distributed_mesh/cell_geometry.F90
           distributed_mesh/cell_topology.F90
           distributed_mesh/dist_mesh_factory.F90
           distributed_mesh/dist_mesh_gmv.F90
	   distributed_mesh/dist_mesh_type.F90
           distributed_mesh/distributed_tet_mesh.F90
           distributed_mesh/facet_hash_type.F90
           distributed_mesh/facet_table_type.F90
           distributed_mesh/unstr_mesh_tools.F90
           distributed_mesh/index_partitioning.F90
           distributed_mesh/mesh_broker.F90
           distributed_mesh/mesh_importer.F90
           distributed_mesh/mesh_modification.F90
           distributed_mesh/parallel_communication.F90
           distributed_mesh/parallel_permutations.F90
           distributed_mesh/simplicial_mesh_support.F90)

set(DISMESH_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/distributed_mesh
	-I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(DISMESH_SOURCE_FILES
                         FILES ${DISMESH_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${DISMESH_FPP_FLAGS}
			 PROCESS_TARGET ${DISMESH_TARGET_NAME})

# Define compile flags
set(DISMESH_COMPILE_FLAGS -I${PGSLib_MODULE_DIR} -I${NETCDF_INCLUDE_DIR})
if(Fortran_COMPILER_IS_INTEL)
  list(APPEND DISMESH_COMPILE_FLAGS "-assume realloc_lhs")
endif()
include(BuildWhitespaceString)
build_whitespace_string(DISMESH_COMPILE_FLAGS_STR ${DISMESH_COMPILE_FLAGS})
set_source_files_properties(${DISMESH_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${DISMESH_COMPILE_FLAGS_STR})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISMESH_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${DISMESH_TARGET_NAME})
	 

