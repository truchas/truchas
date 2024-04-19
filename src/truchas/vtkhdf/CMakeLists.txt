target_sources(truchas PRIVATE
  hdf5_c_binding.F90
  hdf5_c_binding_ext.c
  hl_hdf5.F90
  vtkhdf_file_type.F90
)

add_subdirectory(test)
