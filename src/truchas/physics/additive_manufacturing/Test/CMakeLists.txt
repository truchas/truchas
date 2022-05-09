include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# test_dwell_xyz_motion_type
add_executable(test_toolhead_type test_toolhead_type.F90)
target_link_libraries(test_toolhead_type truchas)

add_test(toolhead_type test_toolhead_type)

set_tests_properties(toolhead_type PROPERTIES PROCESSORS 1)
set_property(TEST toolhead_type PROPERTY LABELS UNIT)
