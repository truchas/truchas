# Truchas files in directory
#   partitioning

# List of files to  process
set(PART_FILES)

# List of files to add to the Truchas library
set(PART_SOURCE_FILES)

# Process target name
set(PART_TARGET_NAME ProcessTruchasPartitioningFiles)

set(PART_FILES
        partitioning/simple_partitioning_methods.F90
        partitioning/graph_partitioner_class.F90
        partitioning/graph_partitioner_factory.F90
        partitioning/chaco_c_binding.F90
        partitioning/chaco_partitioner_type.F90)

set(PART_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        ${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(PART_SOURCE_FILES
                         FILES ${PART_FILES}
                         FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
                         FPP_FLAGS ${PART_FPP_FLAGS}
                         PROCESS_TARGET ${PART_TARGET_NAME})

# Define compile flags
include(BuildWhitespaceString)
set(part_source_flags -I${PGSLib_MODULE_DIR})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND part_source_flags "-assume realloc_lhs")
endif()

build_whitespace_string(part_source_flags_str ${part_source_flags})
set_source_files_properties(${PART_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${part_source_flags_str})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${PART_SOURCE_FILES})
list(APPEND Truchas_PROCESS_TARGETS      ${PART_TARGET_NAME})


