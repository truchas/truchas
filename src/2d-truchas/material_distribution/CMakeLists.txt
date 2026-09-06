target_sources(truchas-extra PRIVATE
  material_distribution_type.F90
  t2d_region_class.F90
  t2d_region_factory.F90
  t2d_cell_set_region_type.F90
  t2d_box_region_type.F90
  t2d_half_plane_region_type.F90
  t2d_disk_region_type.F90
  t2d_region_func_type.F90
  t2d_vol_frac_init_procs.F90
)

add_subdirectory(test)
