include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# test_scalar_func_tools
add_executable(test_scalar_func_tools test_scalar_func_tools.F90)
target_link_libraries(test_scalar_func_tools truchas)
add_test(scalar_func_tools test_scalar_func_tools)
set_tests_properties(scalar_func_tools PROPERTIES PROCESSORS 1)
set_property(TEST scalar_func_tools PROPERTY LABELS UNIT)

if(ENABLE_DYNAMIC_LOADING)
  add_library(test_vector_func_library SHARED test_vector_func_library.F90)

  add_executable(test_vector_func_factory test_vector_func_factory.F90)
  add_dependencies(test_vector_func_factory test_vector_func_library)
  target_link_libraries(test_vector_func_factory truchas)
  target_compile_definitions(test_vector_func_factory PRIVATE
    TEST_VECTOR_FUNC_LIBRARY="$<TARGET_FILE:test_vector_func_library>")
  add_test(vector_func_factory test_vector_func_factory)
  set_tests_properties(vector_func_factory PROPERTIES PROCESSORS 1)
  set_property(TEST vector_func_factory PROPERTY LABELS UNIT)
endif()
