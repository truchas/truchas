# Truchas files in directory
#   grid_mapping

set(GRIDMAP_SOURCE_FILES
           grid_mapping/gm_mesh_type.F90
           grid_mapping/grid_mapping_module.F90
           grid_mapping/grid_mapping_utils.F90
           grid_mapping/hpsort.F90
           grid_mapping/overlap_module.F90)

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${GRIDMAP_SOURCE_FILES})
