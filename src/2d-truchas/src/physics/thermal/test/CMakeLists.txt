include(CTest)

add_library(test_t2d_thermal_common STATIC test_t2d_thermal_common.F90)
set_target_properties(test_t2d_thermal_common PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_include_directories(test_t2d_thermal_common INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(test_t2d_thermal_common PRIVATE ht2d truchas ${Truchas_LINK_LIBS})

set(thermal_unit_tests
  test_t2d_thermal_vector_type
  test_t2d_thermal_model_type
  test_t2d_thermal_precon_type
  test_t2d_thermal_norm_type
  test_t2d_thermal_solver_initial_state)

foreach(target IN LISTS thermal_unit_tests)
  add_executable(${target} ${target}.F90)
  target_link_libraries(${target}
    test_t2d_thermal_common ht2d truchas ${Truchas_LINK_LIBS})
  foreach(nproc IN ITEMS 1 4)
    set(name "${target}-${nproc}p")
    add_test(NAME ${name}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
        $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${name} PROPERTIES
      PROCESSORS ${nproc}
      LABELS "UNIT;CONDUCTION;2d-truchas")
  endforeach()
endforeach()
