include(CTest)

set(flow_unit_tests
  test_t2d_flow_operators
  test_t2d_flow_projection
  test_t2d_flow_projection_solver
  test_t2d_flow_projection_update
  test_t2d_flow_bc_factory
  test_t2d_flow_bc
  test_t2d_flow_momentum
  test_t2d_flow_material_transport
  test_t2d_flow_material_mapping
  test_t2d_flow_material_props
  test_t2d_flow_momentum_solver
  test_t2d_flow_solver
  test_t2d_flow_ic_solver)

foreach(target IN LISTS flow_unit_tests)
  add_executable(${target} ${target}.F90)
  target_link_libraries(${target} flow2d truchas ${Truchas_LINK_LIBS})
  add_test(NAME ${target}-1p
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${target}-1p PROPERTIES
    PROCESSORS 1 LABELS "UNIT;FLOW;2d-truchas")
endforeach()

foreach(target IN ITEMS
    test_t2d_flow_operators
    test_t2d_flow_projection
    test_t2d_flow_projection_solver
    test_t2d_flow_projection_update
    test_t2d_flow_bc_factory
    test_t2d_flow_bc
    test_t2d_flow_momentum
    test_t2d_flow_material_transport
    test_t2d_flow_momentum_solver
    test_t2d_flow_solver
    test_t2d_flow_ic_solver)
  add_test(NAME ${target}-4p
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${target}-4p PROPERTIES
    PROCESSORS 4 LABELS "UNIT;FLOW;2d-truchas")
endforeach()

add_executable(test_t2d_unstr_mesh_factory test_t2d_unstr_mesh_factory.F90)
target_link_libraries(test_t2d_unstr_mesh_factory flow2d truchas ${Truchas_LINK_LIBS})

add_test(NAME flow_operators_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_t2d_flow_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_operators_external-1p PROPERTIES
  PROCESSORS 1 LABELS "UNIT;FLOW;2d-truchas")
add_test(NAME flow_operators_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_t2d_flow_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_operators_external-4p PROPERTIES
  PROCESSORS 4 LABELS "UNIT;FLOW;2d-truchas")

add_test(NAME flow_mesh_factory_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_t2d_unstr_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/mixed-2d.gen
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_mesh_factory_external-1p PROPERTIES
  PROCESSORS 1 LABELS "UNIT;FLOW;2d-truchas")
add_test(NAME flow_mesh_factory_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_t2d_unstr_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_mesh_factory_external-4p PROPERTIES
  PROCESSORS 4 LABELS "UNIT;FLOW;2d-truchas")
