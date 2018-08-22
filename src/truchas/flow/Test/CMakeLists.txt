project(FlowTest Fortran)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/flow)

add_executable(test_scalar_func_copy test_scalar_func_copy.F90)
target_link_libraries(test_scalar_func_copy truchas ${Truchas_LINK_LIBS})

# 1-processor hypre pcg test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${FlowTest_BINARY_DIR}/test_scalar_func_copy ${FlowTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_scalar_func_copy-1p ${test_launch})
set_tests_properties(test_scalar_func_copy-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_scalar_func_copy-1p PROPERTY LABELS UNIT)


# add_executable(test_pressure_poisson test_pressure_poisson.F90)
# target_link_libraries(test_pressure_poisson truchas ${Truchas_LINK_LIBS})

# # 1-processor hypre pcg test
# set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${FlowTest_BINARY_DIR}/test_pressure_poisson ${FlowTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
# add_test(test_pressure_poisson-1p ${test_launch})
# set_tests_properties(test_pressure_poisson-1p PROPERTIES PROCESSORS 1)
# set_property(TEST test_pressure_poisson-1p PROPERTY LABELS UNIT)


add_executable(test_flow_operators test_flow_operators.F90)
target_link_libraries(test_flow_operators truchas ${Truchas_LINK_LIBS})

# 1-processor hypre pcg test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${FlowTest_BINARY_DIR}/test_flow_operators ${FlowTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_flow_operators-1p ${test_launch})
set_tests_properties(test_flow_operators-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_flow_operators-1p PROPERTY LABELS UNIT)
