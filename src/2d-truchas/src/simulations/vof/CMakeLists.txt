add_library(vof2d
  t2d_vof_solver_type.F90
  t2d_vof_sim_type.F90
  t2d_vof_vtkhdf_writer_type.F90)
set_target_properties(vof2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(vof2d PUBLIC volume_tracking_2d truchas-extra truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(vof2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(test)
