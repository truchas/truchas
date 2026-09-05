# --- Library flow2d
set(FLOW2D_SOURCE_FILES
    t2d_flow_state_type.F90
    t2d_flow_operators_type.F90
    t2d_flow_projection_type.F90
    flow_2d_projection_solver_type.F90
    flow_2d_projection_update_type.F90
    t2d_flow_bc_factory_type.F90
    t2d_flow_bc_type.F90
    flow_material_mapping_type.F90
    t2d_flow_material_transport_type.F90
    t2d_flow_momentum_type.F90
    flow_2d_momentum_solver_type.F90
    t2d_flow_material_props_type.F90
    t2d_flow_model_type.F90
    t2d_flow_ic_solver_type.F90
    t2d_flow_solver_type.F90
    t2d_flow_integrator_type.F90
    flow_2d_vtkhdf_writer_type.F90
    t2d_flow_sim_type.F90
)

add_library(flow2d ${FLOW2D_SOURCE_FILES})
set_target_properties(flow2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(flow2d PUBLIC volume_tracking_2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(flow2d BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/../volume_tracking/mod")
target_include_directories(flow2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_subdirectory(Test)
