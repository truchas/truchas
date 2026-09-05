project(HeatTransfer2DTest Fortran)

include(CTest)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)

include_directories(${HeatTransfer2D_BINARY_DIR})

# Compile the shared Fortran fixture once. Adding its source to each test
# executable races to write test_ht_2d_common.mod in a parallel build.
add_library(test_ht_2d_common STATIC support/test_ht_2d_common.F90)
set_target_properties(test_ht_2d_common PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(test_ht_2d_common INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(test_ht_2d_common PRIVATE ht2d truchas ${Truchas_LINK_LIBS})

function(add_ht_unit_test target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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
add_ht_unit_test(test_ht_2d_vector_type test_ht_2d_vector_type.F90)
add_ht_unit_test(test_ht_2d_model_type test_ht_2d_model_type.F90)
add_ht_unit_test(test_ht_2d_precon_type test_ht_2d_precon_type.F90)
add_ht_unit_test(test_ht_2d_norm_type test_ht_2d_norm_type.F90)
add_ht_unit_test(test_ht_2d_solver_initial_state test_ht_2d_solver_initial_state.F90)
add_ht_unit_test(test_ht_2d_solver_type test_ht_2d_solver_type.F90)

function(add_ht_integrated_test name case nproc)
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

function(add_ht_partition_test name case)
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
add_ht_integrated_test(ht_2d_uniform uniform 1)
add_ht_integrated_test(ht_2d_uniform-4p uniform 4)
add_ht_integrated_test(ht_2d_multimaterial_uniform multimaterial_uniform 1)
add_ht_integrated_test(ht_2d_multimaterial_source multimaterial_source 1)
add_ht_integrated_test(ht_2d_multimaterial_mixed_source multimaterial_mixed_source 1)
add_ht_integrated_test(ht_2d_multimaterial_mixed_source-4p multimaterial_mixed_source 4)
add_ht_integrated_test(ht_2d_multimaterial_mixed_tri_source multimaterial_mixed_tri_source 1)
add_ht_integrated_test(ht_2d_multimaterial_conduction multimaterial_conduction 1)
add_ht_integrated_test(ht_2d_linear_quad_x linear 1 test_quad_x.py input_quad_x.json)
add_ht_integrated_test(ht_2d_linear_quad_x-4p linear 4 test_quad_x.py input_quad_x.json)
add_ht_integrated_test(ht_2d_linear_quad_y linear 1 test_quad_y.py input_quad_y.json)
add_ht_integrated_test(ht_2d_linear_quad_y-4p linear 4 test_quad_y.py input_quad_y.json)
add_ht_integrated_test(ht_2d_linear_tri_x linear 1 test_tri_x.py input_tri_x.json)
add_ht_integrated_test(ht_2d_linear_tri_x-4p linear 4 test_tri_x.py input_tri_x.json)
add_ht_integrated_test(ht_2d_linear_tri_y linear 1 test_tri_y.py input_tri_y.json)
add_ht_integrated_test(ht_2d_linear_tri_y-4p linear 4 test_tri_y.py input_tri_y.json)
add_ht_integrated_test(ht_2d_linear_quad_xy linear 1 test_quad_xy.py input_quad_xy.json)
add_ht_integrated_test(ht_2d_linear_quad_xy-4p linear 4 test_quad_xy.py input_quad_xy.json)
add_ht_integrated_test(ht_2d_linear_tri_xy linear 1 test_tri_xy.py input_tri_xy.json)
add_ht_integrated_test(ht_2d_linear_tri_xy-4p linear 4 test_tri_xy.py input_tri_xy.json)
add_ht_integrated_test(ht_2d_linear_mixed_xy linear 1 test_mixed_xy.py input_mixed_xy.json)
add_ht_integrated_test(ht_2d_linear_mixed_xy-4p linear 4 test_mixed_xy.py input_mixed_xy.json)
add_ht_integrated_test(ht_2d_linear_quad_xy_rotated linear 1 test_quad_xy_rotated.py input_quad_xy_rotated.json)
add_ht_integrated_test(ht_2d_linear_quad_xy_rotated-4p linear 4 test_quad_xy_rotated.py input_quad_xy_rotated.json)
add_ht_integrated_test(ht_2d_linear_quad_xy_noisy linear 1 test_quad_xy_noisy.py input_quad_xy_noisy.json)
add_ht_integrated_test(ht_2d_linear_quad_xy_noisy-4p linear 4 test_quad_xy_noisy.py input_quad_xy_noisy.json)
add_ht_integrated_test(ht_2d_source source 1)
add_ht_integrated_test(ht_2d_source-4p source 4)
add_ht_integrated_test(ht_2d_nonlinear nonlinear 1)
add_ht_integrated_test(ht_2d_nonlinear-4p nonlinear 4)
add_ht_integrated_test(ht_2d_nonlinear_inclusion nonlinear_inclusion 4)
add_ht_integrated_test(ht_2d_developer_mesh developer_mesh 4)
add_ht_integrated_test(ht_2d_developer_mesh_tri developer_mesh_tri 4)
add_ht_integrated_test(ht_2d_enthalpy1 enthalpy 1)
add_ht_integrated_test(ht_2d_solidification solidification 4)
add_ht_integrated_test(ht_2d_htc htc 1)
add_ht_integrated_test(ht_2d_htc-4p htc 4)
add_ht_integrated_test(ht_2d_radiation_transient_kelvin radiation 4 test_transient_kelvin.py transient-kelvin.json)
add_ht_integrated_test(ht_2d_radiation_transient_kelvin_scaled radiation 4 test_transient_kelvin_scaled.py transient-kelvin-scaled.json)
add_ht_integrated_test(ht_2d_radiation_transient_celsius radiation 4 test_transient_celsius.py transient-celsius.json)
add_ht_integrated_test(ht_2d_linear_quad_x_htc linear 1 test_quad_x_htc.py input_quad_x_htc.json)
add_ht_integrated_test(ht_2d_linear_quad_x_htc-4p linear 4 test_quad_x_htc.py input_quad_x_htc.json)
add_ht_integrated_test(ht_2d_linear_tri_x_htc linear 1 test_tri_x_htc.py input_tri_x_htc.json)
add_ht_integrated_test(ht_2d_linear_tri_x_htc-4p linear 4 test_tri_x_htc.py input_tri_x_htc.json)

add_ht_partition_test(ht_2d_vtkhdf_partition partition)
add_ht_partition_test(ht_2d_partition_multimaterial_mixed_source multimaterial_mixed_source)
add_ht_partition_test(ht_2d_partition_multimaterial_mixed_tri_source multimaterial_mixed_tri_source)
