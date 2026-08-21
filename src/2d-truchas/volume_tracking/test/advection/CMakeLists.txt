add_vof_test_executable(vof_advection vof_advection.F90)

add_vof_test(circle_advection ${CMAKE_CURRENT_SOURCE_DIR}/test.py vof_advection
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/gold.vtkhdf)
