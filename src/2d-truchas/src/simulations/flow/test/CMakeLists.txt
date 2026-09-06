include(CTest)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)

add_executable(test_t2d_flow_integrator test_t2d_flow_integrator.F90)
target_link_libraries(test_t2d_flow_integrator flow_sim2d truchas ${Truchas_LINK_LIBS})
add_executable(test_t2d_flow_vtkhdf_nan test_t2d_flow_vtkhdf_nan.F90)
target_link_libraries(test_t2d_flow_vtkhdf_nan flow_sim2d truchas ${Truchas_LINK_LIBS})

foreach(nproc IN ITEMS 1 4)
  foreach(target IN ITEMS test_t2d_flow_integrator)
    set(name "${target}-${nproc}p")
    add_test(NAME ${name}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
        $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${name} PROPERTIES
      PROCESSORS ${nproc} LABELS "UNIT;FLOW;2d-truchas")
  endforeach()
endforeach()

add_test(NAME t2d_flow_vtkhdf_nan
  COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_nan/test.py
    $<TARGET_FILE:test_t2d_flow_vtkhdf_nan> ${MPIEXEC})
set_tests_properties(t2d_flow_vtkhdf_nan PROPERTIES
  PROCESSORS 1 TIMEOUT 30 LABELS "UNIT;FLOW;VTKHDF;2d-truchas")

function(add_flow_integrated_test name script input)
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${script}
      $<TARGET_FILE:truchas_2d> ${CMAKE_CURRENT_SOURCE_DIR}/${input}
      ${MPIEXEC} ${ARGN})
  set_tests_properties(${name} PROPERTIES
    PROCESSORS 4 TIMEOUT 120 LABELS "INTEGRATION;FLOW;2d-truchas")
endfunction()

add_flow_integrated_test(flow_stokes_channel_rotated
  stokes/test_flow_integrated.py stokes/pressure_channel_analytic.json)
add_flow_integrated_test(flow_stokes_channel_noisy
  stokes/test_flow_integrated.py stokes/pressure_channel_noisy_analytic.json)
add_flow_integrated_test(flow_channel_rotated
  ns/test_flow_channel.py ns/pressure_channel_analytic.json)
add_flow_integrated_test(flow_stokes_channel
  stokes/test_flow_integrated.py stokes/pressure_channel_analytic_aligned.json 0.0)
add_flow_integrated_test(flow_stokes_couette_viscosity_contrast
  stokes/test_flow_couette.py stokes/couette_viscosity_contrast.json)
add_flow_integrated_test(flow_channel
  ns/test_flow_channel.py ns/pressure_channel_analytic_aligned.json 0.0)
add_flow_integrated_test(flow_inviscid_solid_wall_channel
  ns/test_flow_solid_wall_channel.py ns/solid_wall_channel.json)
add_flow_integrated_test(flow_inviscid_mixed_solid_wall_channel
  ns/test_flow_mixed_solid_wall_channel.py ns/mixed_solid_wall_channel.json)
add_flow_integrated_test(flow_inviscid_solid_wall_hydrostatic_pure
  ns/test_flow_solid_wall_hydrostatic_pure.py ns/solid_wall_hydrostatic_pure.json)
add_flow_integrated_test(flow_inviscid_solid_wall_hydrostatic_mixed
  ns/test_flow_solid_wall_hydrostatic_mixed.py ns/solid_wall_hydrostatic_mixed.json)
add_flow_integrated_test(flow_inviscid_solid_wall_hydrostatic_sloped
  ns/test_flow_solid_wall_hydrostatic_sloped.py ns/solid_wall_hydrostatic_sloped.json)
add_flow_integrated_test(flow_inviscid_solid_wall_hydrostatic_sloped_noisy
  ns/test_flow_solid_wall_hydrostatic_sloped.py ns/solid_wall_hydrostatic_sloped_noisy.json)

# TODO: Enable after the solid/liquid no-slip discretization places the
# effective wall at the material interface rather than at the solid-cell
# center. The current result is not accurate enough for this analytic test.
# add_flow_integrated_test(flow_poiseuille_solid_wall
#   ns/test_flow_poiseuille_solid_wall.py ns/poiseuille_solid_wall.json)

add_flow_integrated_test(flow_couette ns/test_flow_couette.py ns/couette_mixed.json)
add_flow_integrated_test(flow_couette_rotated ns/test_flow_couette.py
  ns/couette_rotated.json 30.0)
add_flow_integrated_test(flow_multifluid_hydrostatic
  ns/test_flow_multifluid_hydrostatic.py ns/multifluid_hydrostatic.json)
add_flow_integrated_test(flow_cylinder ns/cylinder/test.py ns/cylinder/input.json)
add_flow_integrated_test(flow_bubble_rise ns/bubble_rise/test.py ns/bubble_rise/input.json)
add_flow_integrated_test(flow_stokes_body_acceleration
  stokes/test_flow_gravity.py stokes/body_acceleration.json)
add_flow_integrated_test(flow_stokes_hydrostatic
  stokes/test_flow_hydrostatic.py stokes/hydrostatic_box.json)
add_flow_integrated_test(flow_lid_driven_cavity
  ns/test_flow_integrated.py ns/lid_driven_cavity_ns.json)
add_flow_integrated_test(flow_inviscid_channel
  ns/test_flow_inviscid.py ns/inviscid_channel.json)
add_flow_integrated_test(flow_inviscid_channel_rotated
  ns/test_flow_inviscid.py ns/inviscid_channel_rotated.json
  0.8660254037844386 0.5)
add_flow_integrated_test(flow_free_surface_plug
  free_surface/test_pipe.py free_surface/plug.json plug)
add_flow_integrated_test(flow_free_surface_fill
  free_surface/test_pipe.py free_surface/fill.json fill)
add_flow_integrated_test(flow_free_surface_drain
  free_surface/test_pipe.py free_surface/drain.json drain)
