add_library(flow_sim2d
  t2d_flow_integrator_type.F90
  t2d_flow_sim_type.F90
  t2d_flow_vtkhdf_writer_type.F90)
set_target_properties(flow_sim2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(flow_sim2d PUBLIC flow2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(flow_sim2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(test)
