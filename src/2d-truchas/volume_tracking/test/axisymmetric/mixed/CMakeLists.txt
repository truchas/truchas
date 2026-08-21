add_vof_legacy_test(circle_axisymm_mixedmesh vof_axisymmetric
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/reference.txt)
add_vof_vtkhdf_test(circle_axisymm_mixedmesh ${CMAKE_CURRENT_SOURCE_DIR}/test.py vof_axisymmetric
  ${CMAKE_CURRENT_SOURCE_DIR}/input.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/reference.txt)
