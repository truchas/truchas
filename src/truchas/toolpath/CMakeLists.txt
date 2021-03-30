target_sources(truchas PRIVATE
  toolpath_type.F90
  xyz_motion_class.F90
  dwell_xyz_motion_type.F90
  linear_xyz_motion_type.F90
  toolpath_factory.F90
  toolpath_namelist.F90
  toolpath_table.F90
)

add_subdirectory(Test)
