# Truchas files in directory
#   nka

set(NKA_SOURCE_FILES nka/nka_type.F90)

set_source_files_properties(${NKA_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS "-I${PGSLib_MODULE_DIR} -I${Truchas_utilities_dir}")


list(APPEND Truchas_LIBRARY_SOURCE_FILES ${NKA_SOURCE_FILES})
