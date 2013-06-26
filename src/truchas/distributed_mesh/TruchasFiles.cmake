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
           distributed_mesh/distributed_hex_mesh.F90
	   distributed_mesh/distributed_mesh.F90
           distributed_mesh/distributed_mesh_gmv.F90
           distributed_mesh/distributed_tet_mesh.F90
           distributed_mesh/facet_labeling.F90
           distributed_mesh/hashing.F90
           distributed_mesh/hexahedral_mesh_support.F90
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
set_source_files_properties(${DISMESH_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS "-I${PGSLib_MODULE_DIR} -I${NetCDF_INCLUDE_DIR}")

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DISMESH_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${DISMESH_TARGET_NAME})
	 

