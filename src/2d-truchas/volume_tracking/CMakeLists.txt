# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

set(VOLUME_TRACKING_2D_SOURCE_FILES
  geom_axisymmetric.F90
  volume_tracker_2d_class.F90
  simple_volume_tracker_type.F90
  gradient_2d_cc_function.F90
  cell_geom_2d_vof_type.F90
  plane_2d_type.F90
  locate_plane_os_2d_function.F90
  truncation_volume_2d_type.F90
  geometric_volume_tracker_type.F90
  vof_2d_vtkhdf_writer_type.F90
  vof_2d_solver_type.F90
  vof_2d_sim_type.F90
)

add_library(volume_tracking_2d ${VOLUME_TRACKING_2D_SOURCE_FILES})
set_target_properties(volume_tracking_2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(volume_tracking_2d truchas-extra truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(volume_tracking_2d BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(volume_tracking_2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_executable(vof_2d_main vof_2d_main.F90)
target_link_libraries(vof_2d_main volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
set_target_properties(vof_2d_main PROPERTIES OUTPUT_NAME vof_2d)

set(VOLUME_TRACKING_2D_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")
add_subdirectory(test)
