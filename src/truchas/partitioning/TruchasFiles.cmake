# Truchas files in directory
#   partitioning

set(PART_SOURCE_FILES
        partitioning/simple_partitioning_methods.F90
        partitioning/graph_partitioner_class.F90
        partitioning/graph_partitioner_factory.F90
        partitioning/chaco_c_binding.F90
        partitioning/chaco_partitioner_type.F90)

# Define compile flags
include(BuildWhitespaceString)
set(part_source_flags
  -I${PGSLib_MODULE_DIR}
  -I${Truchas_utilities_dir})

build_whitespace_string(part_source_flags_str ${part_source_flags})
set_source_files_properties(${PART_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${part_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${PART_SOURCE_FILES})
