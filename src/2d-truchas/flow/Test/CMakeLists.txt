include(CTest)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)

add_executable(test_t2d_flow_operators test_t2d_flow_operators.F90)
target_link_libraries(test_t2d_flow_operators flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_unstr_mesh_factory test_t2d_unstr_mesh_factory.F90)
target_link_libraries(test_t2d_unstr_mesh_factory flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_projection test_t2d_flow_projection.F90)
target_link_libraries(test_t2d_flow_projection flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_projection_solver test_t2d_flow_projection_solver.F90)
target_link_libraries(test_t2d_flow_projection_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_projection_update test_t2d_flow_projection_update.F90)
target_link_libraries(test_t2d_flow_projection_update flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_bc_factory test_t2d_flow_bc_factory.F90)
target_link_libraries(test_t2d_flow_bc_factory flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_bc test_t2d_flow_bc.F90)
target_link_libraries(test_t2d_flow_bc flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_momentum test_t2d_flow_momentum.F90)
target_link_libraries(test_t2d_flow_momentum flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_material_transport test_t2d_flow_material_transport.F90)
target_link_libraries(test_t2d_flow_material_transport flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_flow_material_mapping test_flow_material_mapping.F90)
target_link_libraries(test_flow_material_mapping flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_material_props test_t2d_flow_material_props.F90)
target_link_libraries(test_t2d_flow_material_props flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_momentum_solver test_t2d_flow_momentum_solver.F90)
target_link_libraries(test_t2d_flow_momentum_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_solver test_t2d_flow_solver.F90)
target_link_libraries(test_t2d_flow_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_integrator test_t2d_flow_integrator.F90)
target_link_libraries(test_t2d_flow_integrator flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_ic_solver test_t2d_flow_ic_solver.F90)
target_link_libraries(test_t2d_flow_ic_solver flow2d truchas ${Truchas_LINK_LIBS})

add_executable(test_t2d_flow_vtkhdf_nan test_t2d_flow_vtkhdf_nan.F90)
target_link_libraries(test_t2d_flow_vtkhdf_nan flow2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_operators ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_operators-1p ${test_launch})
set_tests_properties(t2d_flow_operators-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_operators-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_projection ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection-1p ${test_launch})
set_tests_properties(t2d_flow_projection-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_projection-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_projection_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection_solver-1p ${test_launch})
set_tests_properties(t2d_flow_projection_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_projection_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_projection_update ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection_update-1p ${test_launch})
set_tests_properties(t2d_flow_projection_update-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_projection_update-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_bc_factory ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_bc_factory-1p ${test_launch})
set_tests_properties(t2d_flow_bc_factory-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_bc_factory-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_bc ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_bc-1p ${test_launch})
set_tests_properties(t2d_flow_bc-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_bc-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_momentum ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_momentum-1p ${test_launch})
set_tests_properties(t2d_flow_momentum-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_momentum-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_material_transport ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_material_transport-1p ${test_launch})
set_tests_properties(t2d_flow_material_transport-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_material_transport-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_material_mapping ${MPIEXEC_POSTFLAGS})
add_test(flow_material_mapping-1p ${test_launch})
set_tests_properties(flow_material_mapping-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_material_mapping-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_material_props ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_material_props-1p ${test_launch})
set_tests_properties(t2d_flow_material_props-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_material_props-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_momentum_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_momentum_solver-1p ${test_launch})
set_tests_properties(t2d_flow_momentum_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_momentum_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_solver-1p ${test_launch})
set_tests_properties(t2d_flow_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_integrator ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_integrator-1p ${test_launch})
set_tests_properties(t2d_flow_integrator-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_integrator-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_t2d_flow_ic_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_ic_solver-1p ${test_launch})
set_tests_properties(t2d_flow_ic_solver-1p PROPERTIES PROCESSORS 1)
set_property(TEST t2d_flow_ic_solver-1p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME t2d_flow_vtkhdf_nan
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_nan/test.py
          $<TARGET_FILE:test_t2d_flow_vtkhdf_nan>
          ${MPIEXEC})
set_tests_properties(t2d_flow_vtkhdf_nan PROPERTIES
  PROCESSORS 1
  TIMEOUT 30
  LABELS "UNIT;FLOW;VTKHDF;2d-truchas")

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_operators ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_operators-4p ${test_launch})
set_tests_properties(t2d_flow_operators-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_operators-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME flow_2d_operators_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_t2d_flow_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_operators_external-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_operators_external-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
add_test(NAME flow_2d_operators_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_t2d_flow_operators> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_operators_external-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_operators_external-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME flow_2d_mesh_factory_external-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_t2d_unstr_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/mixed-2d.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_mesh_factory_external-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_mesh_factory_external-1p PROPERTY LABELS UNIT FLOW 2d-truchas)
add_test(NAME flow_2d_mesh_factory_external-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_t2d_unstr_mesh_factory> ${CMAKE_CURRENT_SOURCE_DIR}/unstr-quad.gen ${MPIEXEC_POSTFLAGS})
set_tests_properties(flow_2d_mesh_factory_external-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_mesh_factory_external-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_projection ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection-4p ${test_launch})
set_tests_properties(t2d_flow_projection-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_projection-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_projection_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection_solver-4p ${test_launch})
set_tests_properties(t2d_flow_projection_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_projection_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_projection_update ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_projection_update-4p ${test_launch})
set_tests_properties(t2d_flow_projection_update-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_projection_update-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_bc_factory ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_bc_factory-4p ${test_launch})
set_tests_properties(t2d_flow_bc_factory-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_bc_factory-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_bc ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_bc-4p ${test_launch})
set_tests_properties(t2d_flow_bc-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_bc-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_momentum ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_momentum-4p ${test_launch})
set_tests_properties(t2d_flow_momentum-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_momentum-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_material_transport ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_material_transport-4p ${test_launch})
set_tests_properties(t2d_flow_material_transport-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_material_transport-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_momentum_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_momentum_solver-4p ${test_launch})
set_tests_properties(t2d_flow_momentum_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_momentum_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_solver-4p ${test_launch})
set_tests_properties(t2d_flow_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_integrator ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_integrator-4p ${test_launch})
set_tests_properties(t2d_flow_integrator-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_integrator-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_t2d_flow_ic_solver ${MPIEXEC_POSTFLAGS})
add_test(t2d_flow_ic_solver-4p ${test_launch})
set_tests_properties(t2d_flow_ic_solver-4p PROPERTIES PROCESSORS 4)
set_property(TEST t2d_flow_ic_solver-4p PROPERTY LABELS UNIT FLOW 2d-truchas)

add_test(NAME stokes_2d_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_analytic.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME stokes_2d_integrated_noisy_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_noisy_analytic.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_noisy_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_channel.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/pressure_channel_analytic.json
          ${MPIEXEC})
set_tests_properties(ns_2d_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME stokes_2d_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_integrated.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/pressure_channel_analytic_aligned.json
          ${MPIEXEC} 0.0)
set_tests_properties(stokes_2d_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME stokes_2d_couette_viscosity_contrast
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_couette.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/couette_viscosity_contrast.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_couette_viscosity_contrast PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_channel.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/pressure_channel_analytic_aligned.json
          ${MPIEXEC} 0.0)
set_tests_properties(ns_2d_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_solid_wall_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_solid_wall_channel.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/solid_wall_channel.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_solid_wall_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_mixed_solid_wall_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_mixed_solid_wall_channel.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/mixed_solid_wall_channel.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_mixed_solid_wall_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_solid_wall_hydrostatic_pure
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_solid_wall_hydrostatic_pure.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/solid_wall_hydrostatic_pure.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_solid_wall_hydrostatic_pure PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_solid_wall_hydrostatic_mixed
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_solid_wall_hydrostatic_mixed.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/solid_wall_hydrostatic_mixed.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_solid_wall_hydrostatic_mixed PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_solid_wall_hydrostatic_sloped
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_solid_wall_hydrostatic_sloped.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/solid_wall_hydrostatic_sloped.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_solid_wall_hydrostatic_sloped PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_solid_wall_hydrostatic_sloped_noisy
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_solid_wall_hydrostatic_sloped.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/solid_wall_hydrostatic_sloped_noisy.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_solid_wall_hydrostatic_sloped_noisy PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

# TODO: Enable after the solid/liquid no-slip discretization places the
# effective wall at the material interface rather than at the solid-cell
# center.  The current result is not accurate enough for this analytic test.
# add_test(NAME ns_2d_poiseuille_solid_wall
#   COMMAND ${Python3_EXECUTABLE}
#           ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_poiseuille_solid_wall.py
#           $<TARGET_FILE:truchas_2d>
#           ${CMAKE_CURRENT_SOURCE_DIR}/ns/poiseuille_solid_wall.json
#           ${MPIEXEC})
# set_tests_properties(ns_2d_poiseuille_solid_wall PROPERTIES
#   PROCESSORS 4
#   TIMEOUT 120
#   LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_couette
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_couette.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/couette_mixed.json
          ${MPIEXEC})
set_tests_properties(ns_2d_couette PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_couette_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_couette.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/couette_rotated.json
          ${MPIEXEC}
          30.0)
set_tests_properties(ns_2d_couette_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_multifluid_hydrostatic
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_multifluid_hydrostatic.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/multifluid_hydrostatic.json
          ${MPIEXEC})
set_tests_properties(ns_2d_multifluid_hydrostatic PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_cylinder
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/cylinder/test.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/cylinder/input.json
          ${MPIEXEC})
set_tests_properties(ns_2d_cylinder PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_bubble_rise
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/bubble_rise/test.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/bubble_rise/input.json
          ${MPIEXEC})
set_tests_properties(ns_2d_bubble_rise PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME stokes_2d_integrated_body_acceleration
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_gravity.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/body_acceleration.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_body_acceleration PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME stokes_2d_integrated_hydrostatic
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/test_flow_2d_hydrostatic.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/stokes/hydrostatic_box.json
          ${MPIEXEC})
set_tests_properties(stokes_2d_integrated_hydrostatic PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;STOKES;2d-truchas")

add_test(NAME ns_2d_integrated_lid_cavity
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_integrated.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/lid_driven_cavity_ns.json
          ${MPIEXEC})
set_tests_properties(ns_2d_integrated_lid_cavity PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_inviscid.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/inviscid_channel.json
          ${MPIEXEC})
set_tests_properties(ns_2d_inviscid_channel PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")

add_test(NAME ns_2d_inviscid_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/test_ns_2d_inviscid.py
          $<TARGET_FILE:truchas_2d>
          ${CMAKE_CURRENT_SOURCE_DIR}/ns/inviscid_channel_rotated.json
          ${MPIEXEC} 0.8660254037844386 0.5)
set_tests_properties(ns_2d_inviscid_channel_rotated PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATED;NS;2d-truchas")
