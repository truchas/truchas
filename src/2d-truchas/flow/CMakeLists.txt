# --- Library flow2d
set(FLOW2D_SOURCE_FILES
    flow_2d_state_type.F90
    flow_2d_operators_type.F90
    flow_2d_projection_type.F90
    flow_2d_bc_factory_type.F90
    flow_2d_bc_type.F90
    flow_2d_momentum_type.F90
    flow_2d_momentum_solver_type.F90
)

add_library(flow2d ${FLOW2D_SOURCE_FILES})
set_target_properties(flow2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(flow2d truchas ${Truchas_LINK_LIBS})
target_include_directories(flow2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(Test)
