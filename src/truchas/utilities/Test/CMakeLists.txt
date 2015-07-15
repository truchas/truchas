project(UtilitiesTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
#include_directories(${PGSLib_MODULE_DIR})

add_executable(test_integer_set_type test_integer_set_type.F90)
target_link_libraries(test_integer_set_type truchas ${Truchas_LINK_LIBS})

add_test(integer_set_type test_integer_set_type)

set_tests_properties(integer_set_type PROPERTIES PROCESSORS 1)
set_property(TEST integer_set_type PROPERTY LABELS UNIT)

add_executable(test_graph_type test_graph_type.F90)
target_link_libraries(test_graph_type truchas ${Truchas_LINK_LIBS})

add_test(graph_type test_graph_type)

set_tests_properties(graph_type PROPERTIES PROCESSORS 1)
set_property(TEST graph_type PROPERTY LABELS UNIT)

add_executable(test_permutations test_permutations.F90)
target_link_libraries(test_permutations truchas ${Truchas_LINK_LIBS})

add_test(permutations test_permutations)

set_tests_properties(permutations PROPERTIES PROCESSORS 1)
set_property(TEST permutations PROPERTY LABELS UNIT)
