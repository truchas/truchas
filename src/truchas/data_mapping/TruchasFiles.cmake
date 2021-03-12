# Truchas files in directory
#   data_mapping

set(DATAMAP_SOURCE_FILES
    data_mapping/data_mapper_class.F90
    data_mapping/kuprat_mapper_type.F90
    data_mapping/kuprat/gm_mesh_type.F90
    data_mapping/kuprat/grid_mapping_module.F90
    data_mapping/kuprat/grid_mapping_utils.F90
    data_mapping/kuprat/hpsort.F90
    data_mapping/kuprat/overlap_module.F90
)
if(USE_PORTAGE)
  list(APPEND DATAMAP_SOURCE_FILES data_mapping/portage_mapper_type.F90)
endif()

# Define compile flags
include(BuildWhitespaceString)
build_whitespace_string(datamap_source_flags_string
  -I${Truchas_utilities_dir})
set_source_files_properties(${DATAMAP_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${datamap_source_flags_string})

# Add the C++ files
if(USE_PORTAGE)
  list(APPEND DATAMAP_SOURCE_FILES data_mapping/portage/truchas_portage.cc)
endif()

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DATAMAP_SOURCE_FILES})		


