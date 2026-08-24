# --- Library nsht2d
add_library(nsht2d
  ns_ht_2d_enthalpy_advector_type.F90
  ns_ht_2d_solver_type.F90
  ns_ht_2d_vtkhdf_writer_type.F90
  ns_ht_2d_sim_type.F90)
set_target_properties(nsht2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(nsht2d PUBLIC flow2d ht2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(nsht2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_executable(ns_ht_2d_main ns_ht_2d_main.F90)
target_link_libraries(ns_ht_2d_main nsht2d truchas ${Truchas_LINK_LIBS})
set_target_properties(ns_ht_2d_main PROPERTIES OUTPUT_NAME ns_ht_2d)

add_subdirectory(test)
