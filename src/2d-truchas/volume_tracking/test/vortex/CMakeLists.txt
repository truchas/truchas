add_vof_test_executable(vof_vortex vof_vortex.F90)

add_vof_legacy_test(circle_vortex vof_vortex
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/reference.txt)
add_vof_vtkhdf_test(circle_vortex ${CMAKE_CURRENT_SOURCE_DIR}/test.py vof_vortex
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/reference.txt)
