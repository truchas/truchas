include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test_2d_regions test_2d_regions.F90)
target_link_libraries(test_2d_regions truchas-extra)

add_test(2d_regions test_2d_regions)
set_property(TEST 2d_regions PROPERTY LABELS UNIT)

add_executable(test_region_func_type test_region_func_type.F90)
target_link_libraries(test_region_func_type truchas-extra)

add_test(region_func_type test_region_func_type)
set_property(TEST region_func_type PROPERTY LABELS UNIT)

add_executable(test_compute_volume_fractions test_compute_volume_fractions.F90)
target_link_libraries(test_compute_volume_fractions truchas-extra)

add_test(compute_volume_fractions test_compute_volume_fractions)
set_property(TEST compute_volume_fractions PROPERTY LABELS UNIT)
