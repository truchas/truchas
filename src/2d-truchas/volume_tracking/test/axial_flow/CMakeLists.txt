add_vof_test_executable(vof_axialflow vof_axialflow.F90)

add_vof_test(circle_axialfl ${CMAKE_CURRENT_SOURCE_DIR}/test.py vof_axialflow
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/gold.vtkhdf)
