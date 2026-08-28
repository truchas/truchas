include(CTest)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_executable(test_flow_2d_operators test_flow_2d_operators.F90)
target_link_libraries(test_flow_2d_operators flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_mesh_factory test_flow_2d_mesh_factory.F90)
target_link_libraries(test_flow_2d_mesh_factory flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_projection test_flow_2d_projection.F90)
target_link_libraries(test_flow_2d_projection flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_projection_solver test_flow_2d_projection_solver.F90)
target_link_libraries(test_flow_2d_projection_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_projection_update test_flow_2d_projection_update.F90)
target_link_libraries(test_flow_2d_projection_update flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_bc_factory test_flow_2d_bc_factory.F90)
target_link_libraries(test_flow_2d_bc_factory flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_bc test_flow_2d_bc.F90)
target_link_libraries(test_flow_2d_bc flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_momentum test_flow_2d_momentum.F90)
target_link_libraries(test_flow_2d_momentum flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_material_transport test_flow_2d_material_transport.F90)
target_link_libraries(test_flow_2d_material_transport flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_material_layout test_flow_2d_material_layout.F90)
target_link_libraries(test_flow_2d_material_layout flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_momentum_solver test_flow_2d_momentum_solver.F90)
target_link_libraries(test_flow_2d_momentum_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_solver test_flow_2d_solver.F90)
target_link_libraries(test_flow_2d_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_ns_2d_solver test_ns_2d_solver.F90)
target_link_libraries(test_ns_2d_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_2d_ic_solver test_flow_2d_ic_solver.F90)
target_link_libraries(test_flow_2d_ic_solver flow2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_operators ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_operators-1p ${test_launch})
set_tests_properties(flow_2d_operators-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_operators-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_projection ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection-1p ${test_launch})
set_tests_properties(flow_2d_projection-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_projection-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_projection_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection_solver-1p ${test_launch})
set_tests_properties(flow_2d_projection_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_projection_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_projection_update ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection_update-1p ${test_launch})
set_tests_properties(flow_2d_projection_update-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_projection_update-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_bc_factory ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_bc_factory-1p ${test_launch})
set_tests_properties(flow_2d_bc_factory-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_bc_factory-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_bc ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_bc-1p ${test_launch})
set_tests_properties(flow_2d_bc-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_bc-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_momentum ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_momentum-1p ${test_launch})
set_tests_properties(flow_2d_momentum-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_momentum-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_material_transport ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_material_transport-1p ${test_launch})
set_tests_properties(flow_2d_material_transport-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_material_transport-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_material_layout ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_material_layout-1p ${test_launch})
set_tests_properties(flow_2d_material_layout-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_material_layout-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_momentum_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_momentum_solver-1p ${test_launch})
set_tests_properties(flow_2d_momentum_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_momentum_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_solver-1p ${test_launch})
set_tests_properties(flow_2d_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_ns_2d_solver ${MPIEXEC_POSTFLAGS})
add_test(ns_2d_solver-1p ${test_launch})
set_tests_properties(ns_2d_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST ns_2d_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_ic_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_ic_solver-1p ${test_launch})
set_tests_properties(flow_2d_ic_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_ic_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_operators ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_operators-4p ${test_launch})
set_tests_properties(flow_2d_operators-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_operators-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME flow_2d_operators_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_flow_2d_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_operators_external-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_operators_external-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
add_test(NAME flow_2d_operators_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_flow_2d_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_operators_external-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_operators_external-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME flow_2d_mesh_factory_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_flow_2d_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/mixed-2d.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_mesh_factory_external-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_mesh_factory_external-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
add_test(NAME flow_2d_mesh_factory_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_flow_2d_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_mesh_factory_external-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_mesh_factory_external-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_projection ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection-4p ${test_launch})
set_tests_properties(flow_2d_projection-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_projection-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_projection_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection_solver-4p ${test_launch})
set_tests_properties(flow_2d_projection_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_projection_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_projection_update ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_projection_update-4p ${test_launch})
set_tests_properties(flow_2d_projection_update-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_projection_update-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_bc_factory ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_bc_factory-4p ${test_launch})
set_tests_properties(flow_2d_bc_factory-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_bc_factory-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_bc ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_bc-4p ${test_launch})
set_tests_properties(flow_2d_bc-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_bc-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_momentum ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_momentum-4p ${test_launch})
set_tests_properties(flow_2d_momentum-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_momentum-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_material_transport ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_material_transport-4p ${test_launch})
set_tests_properties(flow_2d_material_transport-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_material_transport-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_momentum_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_momentum_solver-4p ${test_launch})
set_tests_properties(flow_2d_momentum_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_momentum_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_solver-4p ${test_launch})
set_tests_properties(flow_2d_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_ns_2d_solver ${MPIEXEC_POSTFLAGS})
add_test(ns_2d_solver-4p ${test_launch})
set_tests_properties(ns_2d_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST ns_2d_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_ic_solver ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_ic_solver-4p ${test_launch})
set_tests_properties(flow_2d_ic_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_ic_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME stokes_2d_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:flow_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_analytic.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME stokes_2d_integrated_noisy_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:flow_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_noisy_analytic.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_noisy_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_channel.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/pressure_channel_analytic.json
          ${MPIEXEC})
set_tests_properties(ns_2d_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME stokes_2d_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:flow_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_analytic_aligned.json
          ${MPIEXEC} 0.0)
set_tests_properties(stokes_2d_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_channel.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/pressure_channel_analytic_aligned.json
          ${MPIEXEC} 0.0)
set_tests_properties(ns_2d_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_couette
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_couette.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/couette_mixed.json
          ${MPIEXEC})
set_tests_properties(ns_2d_couette PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME stokes_2d_integrated_body_acceleration
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_gravity.py
          $<TARGET_FILE:flow_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/body_acceleration.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_body_acceleration PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME stokes_2d_integrated_hydrostatic
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_hydrostatic.py
          $<TARGET_FILE:flow_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/hydrostatic_box.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_hydrostatic PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_integrated_lid_cavity
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_integrated.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/lid_driven_cavity_ns.json
          ${MPIEXEC})
set_tests_properties(ns_2d_integrated_lid_cavity PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_inviscid.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/inviscid_channel.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_inviscid.py
          $<TARGET_FILE:ns_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/inviscid_channel_rotated.json
          ${MPIEXEC} 0.8660254037844386 0.5)
set_tests_properties(ns_2d_inviscid_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")
