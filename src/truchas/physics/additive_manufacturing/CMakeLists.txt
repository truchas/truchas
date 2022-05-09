target_sources(truchas PRIVATE
  toolhead_driver.F90
  toolhead_event_type.F90
  toolhead_table_type.F90
  toolhead_namelist.F90
  toolhead_type.F90
  laser_irrad_class.F90
  laser_irrad_factory.F90
  beam_laser_irrad_type.F90
  gauss_laser_irrad_type.F90
)

add_subdirectory(Test)
