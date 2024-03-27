include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(vtkhdf_demo vtkhdf_demo.F90)
target_link_libraries(vtkhdf_demo truchas)

add_test(vtkhdf_file_type vtkhdf_demo)
set_tests_properties(vtkhdf_file_type PROPERTIES PROCESSORS 1)
set_property(TEST vtkhdf_file_type PROPERTY LABELS UNIT)
