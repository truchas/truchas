# Truchas files in directory
#   exodus

set(EXOMESH_SOURCE_FILES
           exodus/exodus_c_binding.F90
           exodus/exodus_file_type.F90
           exodus/exodus_mesh_type.F90
           exodus/exodus_mesh_io.F90
           exodus/exodus_mesh_factory.F90
           exodus/exodus_truchas_hack.F90)


# Update the Truchas library file list and targets
set_source_files_properties(${EXOMESH_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS "-I${Truchas_utilities_dir}")
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${EXOMESH_SOURCE_FILES})
