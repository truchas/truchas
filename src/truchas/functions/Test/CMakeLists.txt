include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# test_scalar_func_tools
add_executable(test_scalar_func_tools test_scalar_func_tools.F90)
target_link_libraries(test_scalar_func_tools truchas)
add_test(scalar_func_tools test_scalar_func_tools)
set_tests_properties(scalar_func_tools PROPERTIES PROCESSORS 1)
set_property(TEST scalar_func_tools PROPERTY LABELS UNIT)
