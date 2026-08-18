include(CTest)

add_executable(test_flow_2d_operators test_flow_2d_operators.F90)
target_link_libraries(test_flow_2d_operators flow2d truchas ${Truchas_LINK_LIBS})

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_flow_2d_operators ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_operators-1p ${test_launch})
set_tests_properties(flow_2d_operators-1p PROPERTIES PROCESSORS 1)
set_property(TEST flow_2d_operators-1p PROPERTY LABELS UNIT FLOW 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_flow_2d_operators ${MPIEXEC_POSTFLAGS})
add_test(flow_2d_operators-4p ${test_launch})
set_tests_properties(flow_2d_operators-4p PROPERTIES PROCESSORS 4)
set_property(TEST flow_2d_operators-4p PROPERTY LABELS UNIT FLOW 2d-truchas)
