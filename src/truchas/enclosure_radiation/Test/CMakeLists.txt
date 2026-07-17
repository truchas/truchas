include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test_rade test_rade.F90)
target_link_libraries(test_rade truchas)

set(test_launch-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./test_rade ${CMAKE_CURRENT_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_rade-1p ${test_launch-1p})
set_tests_properties(test_rade-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_rade-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test_rade ${CMAKE_CURRENT_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_rade-4p ${test_launch-4p})
set_tests_properties(test_rade-4p PROPERTIES PROCESSORS 4)
set_property(TEST test_rade-4p PROPERTY LABELS UNIT RADIATION)
