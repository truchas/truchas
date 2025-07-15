include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test_mfd_polygon_matrix_inv test_mfd_polygon_matrix_inv.F90)
target_link_libraries(test_mfd_polygon_matrix_inv truchas)

add_test(mfd_polygon_matrix_inv ${SERIAL_COMMAND_PREFACE}
         ${CMAKE_CURRENT_BINARY_DIR}/test_mfd_polygon_matrix_inv)

set_tests_properties(mfd_polygon_matrix_inv PROPERTIES PROCESSORS 1)
set_property(TEST mfd_polygon_matrix_inv PROPERTY LABELS UNIT)

add_executable(test_multicomp_lever_type test_multicomp_lever_type.F90)
target_link_libraries(test_multicomp_lever_type truchas)

add_executable(run test.F90 block_linear_solver.F90)
target_link_libraries(run truchas)
