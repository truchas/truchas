project(HyprePCGTypeTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})

add_executable(test_hypre_pcg_type test_hypre_pcg_type.F90)
target_link_libraries(test_hypre_pcg_type truchas ${Truchas_LINK_LIBS})

if(ENABLE_MPI)
  set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${HyprePCGTypeTest_BINARY_DIR}/test_hypre_pcg_type ${MPIEXEC_POSTFLAGS})
else()
  set(test_launch test_hypre_pcg_type)
endif()
add_test(hypre_pcg_type ${test_launch})

set_tests_properties(hypre_pcg_type PROPERTIES PROCESSORS 4)
set_property(TEST hypre_pcg_type PROPERTY LABELS UNIT)

add_executable(test_hypre_hybrid_type test_hypre_hybrid_type.F90)
target_link_libraries(test_hypre_hybrid_type truchas ${Truchas_LINK_LIBS})

if(ENABLE_MPI)
  set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${HyprePCGTypeTest_BINARY_DIR}/test_hypre_hybrid_type ${MPIEXEC_POSTFLAGS})
else()
  set(test_launch test_hypre_hybrid_type)
endif()
add_test(hypre_hybrid_type ${test_launch})

set_tests_properties(hypre_hybrid_type PROPERTIES PROCESSORS 4)
set_property(TEST hypre_hybrid_type PROPERTY LABELS UNIT)


