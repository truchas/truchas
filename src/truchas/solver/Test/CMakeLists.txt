include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test_hypre_pcg_type test_hypre_pcg_type.F90)
target_link_libraries(test_hypre_pcg_type truchas)

# 1-processor hypre pcg test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./test_hypre_pcg_type ${MPIEXEC_POSTFLAGS})
add_test(hypre_pcg_type-1p ${test_launch})
set_tests_properties(hypre_pcg_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST hypre_pcg_type-1p PROPERTY LABELS UNIT)

# 4-processor hypre pcg test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test_hypre_pcg_type ${MPIEXEC_POSTFLAGS})
add_test(hypre_pcg_type-4p ${test_launch})
set_tests_properties(hypre_pcg_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST hypre_pcg_type-4p PROPERTY LABELS UNIT)

add_executable(test_hypre_hybrid_type test_hypre_hybrid_type.F90)
target_link_libraries(test_hypre_hybrid_type truchas)

# 1-processor hypre hybrid test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./test_hypre_hybrid_type ${MPIEXEC_POSTFLAGS})
add_test(hypre_hybrid_type-1p ${test_launch})
set_tests_properties(hypre_hybrid_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST hypre_hybrid_type-1p PROPERTY LABELS UNIT)

# 4-processor hypre hybrid test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test_hypre_hybrid_type ${MPIEXEC_POSTFLAGS})
add_test(hypre_hybrid_type-4p ${test_launch})
set_tests_properties(hypre_hybrid_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST hypre_hybrid_type-4p PROPERTY LABELS UNIT)
