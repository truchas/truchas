add_library(thermal_sim2d
  t2d_thermal_integrator_type.F90
  t2d_thermal_sim_type.F90
  t2d_thermal_vtkhdf_writer_type.F90)
set_target_properties(thermal_sim2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(thermal_sim2d PUBLIC ht2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(thermal_sim2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(test)
