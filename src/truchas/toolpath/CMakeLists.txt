target_sources(truchas PRIVATE
  toolpath_type.F90
  xyz_motion_class.F90
  dwell_xyz_motion_type.F90
  linear_xyz_motion_type.F90
  toolpath_factory.F90
  toolpath_factory_type.F90
  toolpath_namelist.F90
  toolpath_driver.F90
  toolpath_event_type.F90
)

add_subdirectory(Test)
