add_library(nsht2d
  t2d_flow_thermal_enthalpy_advector_type.F90
  t2d_flow_thermal_solver_type.F90
  t2d_flow_thermal_integrator_type.F90
  t2d_flow_thermal_vtkhdf_writer_type.F90
  t2d_flow_thermal_sim_type.F90)
set_target_properties(nsht2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(nsht2d PUBLIC thermal_sim2d flow2d ht2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(nsht2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(test)
