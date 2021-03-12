include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})

# test_scalar_func_tools
add_executable(test_scalar_func_tools test_scalar_func_tools.F90)
target_link_libraries(test_scalar_func_tools truchas ${Truchas_LINK_LIBS})
add_test(scalar_func_tools test_scalar_func_tools)
set_tests_properties(scalar_func_tools PROPERTIES PROCESSORS 1)
set_property(TEST scalar_func_tools PROPERTY LABELS UNIT)
