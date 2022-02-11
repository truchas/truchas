# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(reduce_test reduce_test.F90)
target_link_libraries(reduce_test truchas MPI::MPI_Fortran)
add_test(reduce ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./reduce_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(reduce PROPERTIES PROCESSORS 4)
set_property(TEST reduce PROPERTY LABELS UNIT)

add_executable(test-broadcast broadcast_test.F90)
target_link_libraries(test-broadcast truchas MPI::MPI_Fortran)
add_test(broadcast ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test-broadcast ${MPIEXEC_POSTFLAGS})
set_tests_properties(broadcast PROPERTIES PROCESSORS 4)
set_property(TEST broadcast PROPERTY LABELS UNIT)

add_executable(distribute_test distribute_test.F90)
target_link_libraries(distribute_test truchas MPI::MPI_Fortran)
add_test(distribute ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./distribute_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(distribute PROPERTIES PROCESSORS 4)
set_property(TEST distribute PROPERTY LABELS UNIT)

add_executable(collate_test collate_test.F90)
target_link_libraries(collate_test truchas MPI::MPI_Fortran)
add_test(collate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./collate_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(collate PROPERTIES PROCESSORS 4)
set_property(TEST collate PROPERTY LABELS UNIT)

add_executable(test-maxloc maxloc_test.F90)
target_link_libraries(test-maxloc truchas MPI::MPI_Fortran)
add_test(maxloc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ./test-maxloc ${MPIEXEC_POSTFLAGS})
set_tests_properties(maxloc PROPERTIES PROCESSORS 3)
set_property(TEST maxloc PROPERTY LABELS UNIT)
