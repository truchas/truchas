# Truchas files in directory
#   grid_mapping

# List of files to  process
set(GRIDMAP_FILES)

# List of files to add to the Truchas library
set(GRIDMAP_SOURCE_FILES)

# Process target name
set(GRIDMAP_TARGET_NAME ProcessTruchasGridMapFiles)


set(GRIDMAP_FILES
           grid_mapping/gm_mesh_type.F90
           grid_mapping/grid_mapping_module.F90
           grid_mapping/grid_mapping_utils.F90
           grid_mapping/hpsort.F90
           grid_mapping/overlap_module.F90)

set(GRIDMAP_FPP_FLAGS 
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(GRIDMAP_SOURCE_FILES
                         FILES ${GRIDMAP_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${GRIDMAP_FPP_FLAGS}
			 PROCESS_TARGET ${GRIDMAP_TARGET_NAME})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${GRIDMAP_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${GRIDMAP_TARGET_NAME})


