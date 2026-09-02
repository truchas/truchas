# --- Library flow2d
set(FLOW2D_SOURCE_FILES
    flow_2d_state_type.F90
    flow_2d_operators_type.F90
    flow_2d_projection_type.F90
    flow_2d_projection_solver_type.F90
    flow_2d_projection_update_type.F90
    flow_2d_bc_factory_type.F90
    flow_2d_bc_type.F90
    flow_material_mapping_type.F90
    flow_2d_material_transport_type.F90
    flow_2d_momentum_type.F90
    flow_2d_momentum_solver_type.F90
    flow_2d_material_props_type.F90
    flow_2d_model_type.F90
    flow_2d_ic_solver_type.F90
    flow_2d_solver_type.F90
    stokes_2d_solver_type.F90
    ns_2d_solver_type.F90
    flow_2d_vtkhdf_writer_type.F90
    flow_2d_sim_type.F90
    ns_2d_sim_type.F90
)

add_library(flow2d ${FLOW2D_SOURCE_FILES})
set_target_properties(flow2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(flow2d PUBLIC volume_tracking_2d truchas-extra truchas ${Truchas_LINK_LIBS})
target_include_directories(flow2d BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/../volume_tracking/mod")
target_include_directories(flow2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_executable(flow_2d_main flow_2d_main.F90)
target_link_libraries(flow_2d_main flow2d truchas ${Truchas_LINK_LIBS})
set_target_properties(flow_2d_main PROPERTIES OUTPUT_NAME flow_2d)

add_executable(ns_2d_main ns_2d_main.F90)
target_link_libraries(ns_2d_main flow2d truchas ${Truchas_LINK_LIBS})
set_target_properties(ns_2d_main PROPERTIES OUTPUT_NAME ns_2d)

add_subdirectory(Test)
