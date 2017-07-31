include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})

# test_dwell_xyz_motion_type
add_executable(test_ded_head_type test_ded_head_type.F90)
target_link_libraries(test_ded_head_type truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_ded_head_type.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()

add_test(ded_head_type test_ded_head_type)

set_tests_properties(ded_head_type PROPERTIES PROCESSORS 1)
set_property(TEST ded_head_type PROPERTY LABELS UNIT)
