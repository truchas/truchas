# Truchas files in directory
#   communication

set(COMM_SOURCE_FILES
         communication/parallel_info_module.F90
         communication/parallel_util_module.F90)

set_source_files_properties(${COMM_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS -I${PGSLib_MODULE_DIR})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${COMM_SOURCE_FILES})
