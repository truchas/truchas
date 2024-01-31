target_sources(truchas-extra PRIVATE
  region_class.F90
  region_factory.F90
  cell_set_region_type.F90
  box_region_type.F90
  half_plane_region_type.F90
  disk_region_type.F90
  region_func_type.F90
  vol_frac_init_procs.F90
)

add_subdirectory(test)
