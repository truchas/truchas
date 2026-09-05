add_library(vof_vortex_velocity SHARED vortex_velocity.F90)

file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input.json"
  INPUT "${CMAKE_CURRENT_SOURCE_DIR}/input.json.in")

add_vof_test(circle_vortex ${CMAKE_CURRENT_SOURCE_DIR}/test.py
  ${CMAKE_CURRENT_BINARY_DIR}/input.json
  ${CMAKE_CURRENT_SOURCE_DIR}/gold.vtkhdf)
