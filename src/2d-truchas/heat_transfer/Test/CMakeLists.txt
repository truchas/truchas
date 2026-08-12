project(HeatTransfer2DTest Fortran)

include(CTest)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)

include_directories(${HeatTransfer2D_BINARY_DIR})

# Compile the shared Fortran module once.  Adding its source to each test
# executable races to write test_ht_2d_common.mod in a parallel build.
add_library(test_ht_2d_common STATIC test_HT_2d_common.F90)
set_target_properties(test_ht_2d_common PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(test_ht_2d_common INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(test_ht_2d_common PRIVATE ht2d truchas ${Truchas_LINK_LIBS})


# --- ht_2d_vector_type tests
add_executable(test_ht_2d_vector_type test_ht_2d_vector_type.F90)
target_link_libraries(test_ht_2d_vector_type ht2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_ht_2d_vector_type ${MPIEXEC_POSTFLAGS})
add_test(ht_2d_vector_type-1p ${test_launch})
set_tests_properties(ht_2d_vector_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST ht_2d_vector_type-1p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_ht_2d_vector_type ${MPIEXEC_POSTFLAGS})
add_test(ht_2d_vector_type-4p ${test_launch})
set_tests_properties(ht_2d_vector_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST ht_2d_vector_type-4p PROPERTY LABELS UNIT CONDUCTION 2d-truchas)


# --- HT_2D_MODEL_TYPE tests
add_executable(test_HT_2d_model_type test_HT_2d_model_type.F90)
target_link_libraries(test_HT_2d_model_type test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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
add_executable(test_HT_2d_precon_type test_HT_2d_precon_type.F90)
target_link_libraries(test_HT_2d_precon_type test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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
add_executable(test_HT_2d_norm_type test_HT_2d_norm_type.F90)
target_link_libraries(test_HT_2d_norm_type test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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
add_executable(test_HT_2d_solver_initial_state test_HT_2d_solver_initial_state.F90)
target_link_libraries(test_HT_2d_solver_initial_state test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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
add_executable(test_HT_2d_solver_type test_HT_2d_solver_type.F90)
target_link_libraries(test_HT_2d_solver_type test_ht_2d_common ht2d truchas ${Truchas_LINK_LIBS})

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

add_test(NAME HT_2d_integrated_source
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/source.json
          source 1 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_source PROPERTIES
  PROCESSORS 1
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_nonlinear
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/nonlinear.json
          nonlinear)
set_tests_properties(HT_2d_integrated_nonlinear PROPERTIES
  PROCESSORS 1
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_developer_mesh-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/developer_mesh.json
          developer_mesh 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_developer_mesh-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_developer_mesh_tri-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/developer_mesh_tri.json
          developer_mesh_tri 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_developer_mesh_tri-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_uniform-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/uniform.json
          uniform 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_uniform-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_linear-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/linear.json
          linear 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_linear-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_source-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/source.json
          source 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_source-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")

add_test(NAME HT_2d_integrated_nonlinear-4p
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/test_ht_2d_integrated.py
          $<TARGET_FILE:HT_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/nonlinear.json
          nonlinear 4 ${MPIEXEC})
set_tests_properties(HT_2d_integrated_nonlinear-4p PROPERTIES
  PROCESSORS 4
  TIMEOUT 120
  LABELS "INTEGRATION;CONDUCTION;2d-truchas")
