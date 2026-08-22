add_vof_test(circle_advection ${CMAKE_CURRENT_SOURCE_DIR}/test.py vof_2d_main
  ${CMAKE_CURRENT_SOURCE_DIR}/input.json
  ${CMAKE_CURRENT_SOURCE_DIR}/gold.vtkhdf)
