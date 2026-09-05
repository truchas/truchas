# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

set(VOLUME_TRACKING_2D_SOURCE_FILES
  t2d_geom_axisymmetric.F90
  t2d_volume_tracker_class.F90
  t2d_simple_volume_tracker_type.F90
  t2d_gradient_cc_function.F90
  t2d_cell_geom_vof_type.F90
  t2d_plane_type.F90
  t2d_locate_plane_os_function.F90
  t2d_truncation_volume_type.F90
  t2d_geometric_volume_tracker_type.F90
  t2d_vof_vtkhdf_writer_type.F90
  t2d_vof_solver_type.F90
  t2d_vof_sim_type.F90
)

add_library(volume_tracking_2d ${VOLUME_TRACKING_2D_SOURCE_FILES})
set_target_properties(volume_tracking_2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(volume_tracking_2d truchas-extra truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(volume_tracking_2d BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(volume_tracking_2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

set(VOLUME_TRACKING_2D_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")
add_subdirectory(test)
