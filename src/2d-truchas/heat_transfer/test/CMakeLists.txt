project(HeatTransfer2DTest Fortran)

include(CTest)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)

include_directories(${HeatTransfer2D_BINARY_DIR})

# Compile the shared Fortran fixture once. Adding its source to each test
# executable races to write test_t2d_thermal_common.mod in a parallel build.
add_library(test_t2d_thermal_common STATIC support/test_t2d_thermal_common.F90)
set_target_properties(test_t2d_thermal_common PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(test_t2d_thermal_common INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(test_t2d_thermal_common PRIVATE ht2d truchas ${Truchas_LINK_LIBS})

function(add_thermal_unit_test target source)
  add_executable(${target} ${source})
  if(target STREQUAL "test_t2d_thermal_integrator")
    target_link_libraries(${target} test_t2d_thermal_common thermal_sim2d truchas ${Truchas_LINK_LIBS})
  else()
    target_link_libraries(${target} test_t2d_thermal_common ht2d truchas ${Truchas_LINK_LIBS})
  endif()

  foreach(nproc IN ITEMS 1 4)
    set(name "${target}-${nproc}p")
    add_test(NAME ${name}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
        $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${name} PROPERTIES
      PROCESSORS ${nproc}
      LABELS "UNIT;CONDUCTION;2d-truchas")
  endforeach()
endfunction()

# Component-level tests.
add_thermal_unit_test(test_t2d_thermal_vector_type test_t2d_thermal_vector_type.F90)
add_thermal_unit_test(test_t2d_thermal_model_type test_t2d_thermal_model_type.F90)
add_thermal_unit_test(test_t2d_thermal_precon_type test_t2d_thermal_precon_type.F90)
add_thermal_unit_test(test_t2d_thermal_norm_type test_t2d_thermal_norm_type.F90)
add_thermal_unit_test(test_t2d_thermal_solver_initial_state test_t2d_thermal_solver_initial_state.F90)
add_thermal_unit_test(test_t2d_thermal_integrator test_t2d_thermal_integrator.F90)

function(add_thermal_integrated_test name case nproc)
  set(script_name test.py)
  set(input_name input.json)
  if(ARGC GREATER 3)
    set(script_name ${ARGV3})
  endif()
  if(ARGC GREATER 4)
    set(input_name ${ARGV4})
  endif()
  set(script ${CMAKE_CURRENT_SOURCE_DIR}/${case}/${script_name})
  set(input ${CMAKE_CURRENT_SOURCE_DIR}/${case}/${input_name})
  if(nproc GREATER 1)
    add_test(NAME ${name}
      COMMAND ${Python3_EXECUTABLE} ${script}
        $<TARGET_FILE:truchas_2d> ${nproc} ${MPIEXEC})
  else()
    add_test(NAME ${name}
      COMMAND ${Python3_EXECUTABLE} ${script}
        $<TARGET_FILE:truchas_2d>)
  endif()
  set_tests_properties(${name} PROPERTIES
    PROCESSORS ${nproc}
    TIMEOUT 120
    LABELS "INTEGRATION;CONDUCTION;2d-truchas")
endfunction()

function(add_thermal_partition_test name case)
  set(script_name partition_test.py)
  if(ARGC GREATER 2)
    set(script_name ${ARGV2})
  endif()
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/${case}/${script_name}
      $<TARGET_FILE:truchas_2d> ${MPIEXEC})
  set_tests_properties(${name} PROPERTIES
    PROCESSORS 4
    TIMEOUT 120
    LABELS "INTEGRATION;CONDUCTION;2d-truchas")
endfunction()

# JSON-driven application tests. Each case owns its input and Python checker.
add_thermal_integrated_test(thermal_uniform uniform 1)
add_thermal_integrated_test(thermal_uniform-4p uniform 4)
add_thermal_integrated_test(thermal_multimaterial_uniform multimaterial_uniform 1)
add_thermal_integrated_test(thermal_multimaterial_source multimaterial_source 1)
add_thermal_integrated_test(thermal_multimaterial_mixed_source multimaterial_mixed_source 1)
add_thermal_integrated_test(thermal_multimaterial_mixed_source-4p multimaterial_mixed_source 4)
add_thermal_integrated_test(thermal_multimaterial_mixed_tri_source multimaterial_mixed_tri_source 1)
add_thermal_integrated_test(thermal_multimaterial_conduction multimaterial_conduction 1)
add_thermal_integrated_test(thermal_linear_quad_x linear 1 test_quad_x.py input_quad_x.json)
add_thermal_integrated_test(thermal_linear_quad_x-4p linear 4 test_quad_x.py input_quad_x.json)
add_thermal_integrated_test(thermal_linear_quad_y linear 1 test_quad_y.py input_quad_y.json)
add_thermal_integrated_test(thermal_linear_quad_y-4p linear 4 test_quad_y.py input_quad_y.json)
add_thermal_integrated_test(thermal_linear_tri_x linear 1 test_tri_x.py input_tri_x.json)
add_thermal_integrated_test(thermal_linear_tri_x-4p linear 4 test_tri_x.py input_tri_x.json)
add_thermal_integrated_test(thermal_linear_tri_y linear 1 test_tri_y.py input_tri_y.json)
add_thermal_integrated_test(thermal_linear_tri_y-4p linear 4 test_tri_y.py input_tri_y.json)
add_thermal_integrated_test(thermal_linear_quad_xy linear 1 test_quad_xy.py input_quad_xy.json)
add_thermal_integrated_test(thermal_linear_quad_xy-4p linear 4 test_quad_xy.py input_quad_xy.json)
add_thermal_integrated_test(thermal_linear_tri_xy linear 1 test_tri_xy.py input_tri_xy.json)
add_thermal_integrated_test(thermal_linear_tri_xy-4p linear 4 test_tri_xy.py input_tri_xy.json)
add_thermal_integrated_test(thermal_linear_mixed_xy linear 1 test_mixed_xy.py input_mixed_xy.json)
add_thermal_integrated_test(thermal_linear_mixed_xy-4p linear 4 test_mixed_xy.py input_mixed_xy.json)
add_thermal_integrated_test(thermal_linear_quad_xy_rotated linear 1 test_quad_xy_rotated.py input_quad_xy_rotated.json)
add_thermal_integrated_test(thermal_linear_quad_xy_rotated-4p linear 4 test_quad_xy_rotated.py input_quad_xy_rotated.json)
add_thermal_integrated_test(thermal_linear_quad_xy_noisy linear 1 test_quad_xy_noisy.py input_quad_xy_noisy.json)
add_thermal_integrated_test(thermal_linear_quad_xy_noisy-4p linear 4 test_quad_xy_noisy.py input_quad_xy_noisy.json)
add_thermal_integrated_test(thermal_source source 1)
add_thermal_integrated_test(thermal_source-4p source 4)
add_thermal_integrated_test(thermal_nonlinear nonlinear 1)
add_thermal_integrated_test(thermal_nonlinear-4p nonlinear 4)
add_thermal_integrated_test(thermal_nonlinear_inclusion nonlinear_inclusion 4)
add_thermal_integrated_test(thermal_developer_mesh developer_mesh 4)
add_thermal_integrated_test(thermal_developer_mesh_tri developer_mesh_tri 4)
add_thermal_integrated_test(thermal_enthalpy1 enthalpy 1)
add_thermal_integrated_test(thermal_solidification solidification 4)
add_thermal_integrated_test(thermal_htc htc 1)
add_thermal_integrated_test(thermal_htc-4p htc 4)
add_thermal_integrated_test(thermal_radiation_transient_kelvin radiation 4 test_transient_kelvin.py transient-kelvin.json)
add_thermal_integrated_test(thermal_radiation_transient_kelvin_scaled radiation 4 test_transient_kelvin_scaled.py transient-kelvin-scaled.json)
add_thermal_integrated_test(thermal_radiation_transient_celsius radiation 4 test_transient_celsius.py transient-celsius.json)
add_thermal_integrated_test(thermal_linear_quad_x_htc linear 1 test_quad_x_htc.py input_quad_x_htc.json)
add_thermal_integrated_test(thermal_linear_quad_x_htc-4p linear 4 test_quad_x_htc.py input_quad_x_htc.json)
add_thermal_integrated_test(thermal_linear_tri_x_htc linear 1 test_tri_x_htc.py input_tri_x_htc.json)
add_thermal_integrated_test(thermal_linear_tri_x_htc-4p linear 4 test_tri_x_htc.py input_tri_x_htc.json)

add_thermal_partition_test(thermal_vtkhdf_partition partition)
add_thermal_partition_test(thermal_partition_multimaterial_mixed_source multimaterial_mixed_source)
add_thermal_partition_test(thermal_partition_multimaterial_mixed_tri_source multimaterial_mixed_tri_source)
