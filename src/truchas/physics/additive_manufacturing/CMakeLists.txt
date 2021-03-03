target_sources(truchas PRIVATE
  ded_head_driver.F90
  ded_head_namelist.F90
  ded_head_type.F90
  laser_irrad_class.F90
  laser_irrad_factory.F90
  beam_laser_irrad_type.F90
  gauss_laser_irrad_type.F90
)

add_subdirectory(Test)
