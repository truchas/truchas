include(CTest)

add_executable(test_simulation_log test_simulation_log.F90)
target_link_libraries(test_simulation_log truchas-extra truchas ${Truchas_LINK_LIBS})

add_test(NAME simulation_log-1p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  $<TARGET_FILE:test_simulation_log> ${MPIEXEC_POSTFLAGS})
set_tests_properties(simulation_log-1p PROPERTIES PROCESSORS 1)
set_property(TEST simulation_log-1p PROPERTY LABELS UNIT 2d-truchas)

add_test(NAME simulation_log-4p COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  $<TARGET_FILE:test_simulation_log> ${MPIEXEC_POSTFLAGS})
set_tests_properties(simulation_log-4p PROPERTIES PROCESSORS 4)
set_property(TEST simulation_log-4p PROPERTY LABELS UNIT 2d-truchas)
