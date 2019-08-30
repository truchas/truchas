project(ToolsRadeTest Fortran)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasToolsRade_BINARY_DIR})

add_executable(test_rade_tools test_rade_tools.F90)
target_link_libraries(test_rade_tools rade truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_rade_tools.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()


set(test_launch-synthetic-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${ToolsRadeTest_BINARY_DIR}/test_rade_tools
  ${MPIEXEC_POSTFLAGS})
add_test(test_rade_tools-synthetic-1p ${test_launch-synthetic-1p})
set_tests_properties(test_rade_tools-synthetic-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_rade_tools-synthetic-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${ToolsRadeTest_BINARY_DIR}/test_rade_tools
  ${ToolsRadeTest_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(test_rade_tools-face-1p ${test_launch-face-1p})
set_tests_properties(test_rade_tools-face-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_rade_tools-face-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-face-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${ToolsRadeTest_BINARY_DIR}/test_rade_tools
  ${ToolsRadeTest_SOURCE_DIR}/test_rade_tools_face.nc ${MPIEXEC_POSTFLAGS})
add_test(test_rade_tools-face-4p ${test_launch-face-4p})
set_tests_properties(test_rade_tools-face-4p PROPERTIES PROCESSORS 4)
set_property(TEST test_rade_tools-face-4p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${ToolsRadeTest_BINARY_DIR}/test_rade_tools
  ${ToolsRadeTest_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(test_rade_tools-patch-1p ${test_launch-patch-1p})
set_tests_properties(test_rade_tools-patch-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_rade_tools-patch-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-patch-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${ToolsRadeTest_BINARY_DIR}/test_rade_tools
  ${ToolsRadeTest_SOURCE_DIR}/test_rade_tools_patch.nc ${MPIEXEC_POSTFLAGS})
add_test(test_rade_tools-patch-4p ${test_launch-patch-4p})
set_tests_properties(test_rade_tools-patch-4p PROPERTIES PROCESSORS 4)
set_property(TEST test_rade_tools-patch-4p PROPERTY LABELS UNIT RADIATION)
