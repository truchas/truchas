include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test_integer_set_type test_integer_set_type.F90)
target_link_libraries(test_integer_set_type truchas)

add_test(integer_set_type ${SERIAL_COMMAND_PREFACE}
         ${CMAKE_CURRENT_BINARY_DIR}/test_integer_set_type)

set_tests_properties(integer_set_type PROPERTIES PROCESSORS 1)
set_property(TEST integer_set_type PROPERTY LABELS UNIT)

add_executable(test_graph_type test_graph_type.F90)
target_link_libraries(test_graph_type truchas)

add_test(graph_type ${SERIAL_COMMAND_PREFACE}
         ${CMAKE_CURRENT_BINARY_DIR}/test_graph_type)

set_tests_properties(graph_type PROPERTIES PROCESSORS 1)
set_property(TEST graph_type PROPERTY LABELS UNIT)

add_executable(test_permutations test_permutations.F90)
target_link_libraries(test_permutations truchas)

add_test(permutations ${SERIAL_COMMAND_PREFACE}
         ${CMAKE_CURRENT_BINARY_DIR}/test_permutations)

set_tests_properties(permutations PROPERTIES PROCESSORS 1)
set_property(TEST permutations PROPERTY LABELS UNIT)
