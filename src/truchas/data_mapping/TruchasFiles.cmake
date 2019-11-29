# Truchas files in directory
#   data_mapping

set(DATAMAP_SOURCE_FILES
    data_mapping/kuprat_mapper_type.F90
    data_mapping/kuprat/gm_mesh_type.F90
    data_mapping/kuprat/grid_mapping_module.F90
    data_mapping/kuprat/grid_mapping_utils.F90
    data_mapping/kuprat/hpsort.F90
    data_mapping/kuprat/overlap_module.F90
)

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DATAMAP_SOURCE_FILES})		       


