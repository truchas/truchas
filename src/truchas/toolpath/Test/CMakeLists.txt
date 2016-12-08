project(ToolpathTest Fortran)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})

# test_dwell_xyz_motion_type
add_executable(test_dwell_xyz_motion_type test_dwell_xyz_motion_type.F90)
target_link_libraries(test_dwell_xyz_motion_type truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_dwell_xyz_motion_type.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()

add_test(dwell_xyz_motion_type test_dwell_xyz_motion_type)

set_tests_properties(dwell_xyz_motion_type PROPERTIES PROCESSORS 1)
set_property(TEST dwell_xyz_motion_type PROPERTY LABELS UNIT)

# test_linear_xyz_motion_type
add_executable(test_linear_xyz_motion_type test_linear_xyz_motion_type.F90)
target_link_libraries(test_linear_xyz_motion_type truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_linear_xyz_motion_type.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()

add_test(linear_xyz_motion_type test_linear_xyz_motion_type)

set_tests_properties(linear_xyz_motion_type PROPERTIES PROCESSORS 1)
set_property(TEST linear_xyz_motion_type PROPERTY LABELS UNIT)

# test_toolpath_type
add_executable(test_toolpath_type test_toolpath_type.F90)
target_link_libraries(test_toolpath_type truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_toolpath_type.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()

add_test(toolpath_type test_toolpath_type)

set_tests_properties(toolpath_type PROPERTIES PROCESSORS 1)
set_property(TEST toolpath_type PROPERTY LABELS UNIT)

# test_toolpath_factory
add_executable(test_toolpath_factory test_toolpath_factory.F90)
target_link_libraries(test_toolpath_factory truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_toolpath_factory.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()

add_test(toolpath_factory test_toolpath_factory ${ToolpathTest_SOURCE_DIR})

set_tests_properties(toolpath_factory PROPERTIES PROCESSORS 1)
set_property(TEST toolpath_factory PROPERTY LABELS UNIT)
