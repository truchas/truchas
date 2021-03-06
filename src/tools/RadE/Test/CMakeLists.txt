include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${TruchasToolsRade_BINARY_DIR})

# --- Common modules
add_library(test_rade_tools OBJECT test_rade_tools_common.F90)

# --- RE_DIST_VF_TYPE tests
add_executable(test_dist_vf_type test_dist_vf_type.F90)
target_link_libraries(test_dist_vf_type test_rade_tools rade truchas)

set(test_launch-synthetic-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_dist_vf_type ${MPIEXEC_POSTFLAGS})
add_test(dist_vf_type-synthetic-1p ${test_launch-synthetic-1p})
set_tests_properties(dist_vf_type-synthetic-1p PROPERTIES PROCESSORS 1)
set_property(TEST dist_vf_type-synthetic-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_dist_vf_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(dist_vf_type-face-1p ${test_launch-face-1p})
set_tests_properties(dist_vf_type-face-1p PROPERTIES PROCESSORS 1)
set_property(TEST dist_vf_type-face-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_dist_vf_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(dist_vf_type-face-4p ${test_launch-face-4p})
set_tests_properties(dist_vf_type-face-4p PROPERTIES PROCESSORS 4)
set_property(TEST dist_vf_type-face-4p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_dist_vf_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(dist_vf_type-patch-1p ${test_launch-patch-1p})
set_tests_properties(dist_vf_type-patch-1p PROPERTIES PROCESSORS 1)
set_property(TEST dist_vf_type-patch-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_dist_vf_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(dist_vf_type-patch-4p ${test_launch-patch-4p})
set_tests_properties(dist_vf_type-patch-4p PROPERTIES PROCESSORS 4)
set_property(TEST dist_vf_type-patch-4p PROPERTY LABELS UNIT RADIATION)


# --- RE_VF_DIFF_TYPE tests
add_executable(test_vf_diff_type test_vf_diff_type.F90)
target_link_libraries(test_vf_diff_type test_rade_tools rade truchas)

set(test_launch-synthetic-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_vf_diff_type
  ${MPIEXEC_POSTFLAGS})
add_test(vf_diff_type-synthetic-1p ${test_launch-synthetic-1p})
set_tests_properties(vf_diff_type-synthetic-1p PROPERTIES PROCESSORS 1)
set_property(TEST vf_diff_type-synthetic-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_vf_diff_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(vf_diff_type-face-1p ${test_launch-face-1p})
set_tests_properties(vf_diff_type-face-1p PROPERTIES PROCESSORS 1)
set_property(TEST vf_diff_type-face-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_vf_diff_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(vf_diff_type-face-4p ${test_launch-face-4p})
set_tests_properties(vf_diff_type-face-4p PROPERTIES PROCESSORS 4)
set_property(TEST vf_diff_type-face-4p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_vf_diff_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(vf_diff_type-patch-1p ${test_launch-patch-1p})
set_tests_properties(vf_diff_type-patch-1p PROPERTIES PROCESSORS 1)
set_property(TEST vf_diff_type-patch-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_vf_diff_type
  ${CMAKE_CURRENT_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(vf_diff_type-patch-4p ${test_launch-patch-4p})
set_tests_properties(vf_diff_type-patch-4p PROPERTIES PROCESSORS 4)
set_property(TEST vf_diff_type-patch-4p PROPERTY LABELS UNIT RADIATION)
