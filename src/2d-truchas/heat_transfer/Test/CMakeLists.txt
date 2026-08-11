project(HeatTransfer2DTest Fortran)

include(CTest)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

include_directories(${HeatTransfer2D_BINARY_DIR})


# --- HT_2D_MODEL_TYPE tests
add_executable(test_HT_2d_model_type test_HT_2d_model_type.F90
                                     test_HT_2d_common.F90)
target_link_libraries(test_HT_2d_model_type ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_HT_2d_model_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_model_type-1p ${test_launch})
set_tests_properties(HT_2d_model_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST HT_2d_model_type-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_HT_2d_model_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_model_type-4p ${test_launch})
set_tests_properties(HT_2d_model_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST HT_2d_model_type-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- HT_2D_PRECON_TYPE tests
add_executable(test_HT_2d_precon_type test_HT_2d_precon_type.F90
                                      test_HT_2d_common.F90)
target_link_libraries(test_HT_2d_precon_type ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_HT_2d_precon_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_precon_type-1p ${test_launch})
set_tests_properties(HT_2d_precon_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST HT_2d_precon_type-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_HT_2d_precon_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_precon_type-4p ${test_launch})
set_tests_properties(HT_2d_precon_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST HT_2d_precon_type-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- HT_2D_NORM_TYPE tests
add_executable(test_HT_2d_norm_type test_HT_2d_norm_type.F90
                                    test_HT_2d_common.F90)
target_link_libraries(test_HT_2d_norm_type ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_HT_2d_norm_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_norm_type-1p ${test_launch})
set_tests_properties(HT_2d_norm_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST HT_2d_norm_type-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_HT_2d_norm_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_norm_type-4p ${test_launch})
set_tests_properties(HT_2d_norm_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST HT_2d_norm_type-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- HT_2D_SOLVER_INITIAL_STATE tests
add_executable(test_HT_2d_solver_initial_state test_HT_2d_solver_initial_state.F90
                                       test_HT_2d_common.F90)
target_link_libraries(test_HT_2d_solver_initial_state ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_HT_2d_solver_initial_state ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_solver_initial_state-1p ${test_launch})
set_tests_properties(HT_2d_solver_initial_state-1p PROPERTIES PROCESSORS 1)
set_property(TEST HT_2d_solver_initial_state-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_HT_2d_solver_initial_state ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_solver_initial_state-4p ${test_launch})
set_tests_properties(HT_2d_solver_initial_state-4p PROPERTIES PROCESSORS 4)
set_property(TEST HT_2d_solver_initial_state-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- HT_2D_SOLVER_TYPE tests
add_executable(test_HT_2d_solver_type test_HT_2d_solver_type.F90
                                      test_HT_2d_common.F90)
target_link_libraries(test_HT_2d_solver_type ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_HT_2d_solver_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_solver_type-1p ${test_launch})
set_tests_properties(HT_2d_solver_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST HT_2d_solver_type-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_HT_2d_solver_type ${MPIEXEC_POSTFLAGS})
add_test(HT_2d_solver_type-4p ${test_launch})
set_tests_properties(HT_2d_solver_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST HT_2d_solver_type-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- End-to-end JSON-driven application tests
add_test(NAME HT_2d_integrated_uniform
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/uniform.json
          uniform)
set_tests_properties(HT_2d_integrated_uniform PROPERTIES
  PROCESSORS 1
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_linear
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/linear.json
          linear)
set_tests_properties(HT_2d_integrated_linear PROPERTIES
  PROCESSORS 1
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")
