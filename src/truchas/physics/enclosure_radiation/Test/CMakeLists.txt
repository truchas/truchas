project(RadETest Fortran)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})

add_executable(test_rade test_rade.F90)
target_link_libraries(test_rade truchas ${Truchas_LINK_LIBS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties(test_rade.F90 PROPERTIES
      COMPILE_FLAGS "-assume realloc_lhs")
endif()


set(test_launch-1p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${RadETest_BINARY_DIR}/test_rade ${RadETest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_rade-1p ${test_launch-1p})
set_tests_properties(test_rade-1p PROPERTIES PROCESSORS 1)
set_property(TEST test_rade-1p PROPERTY LABELS UNIT RADIATION)

set(test_launch-4p ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${RadETest_BINARY_DIR}/test_rade ${RadETest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(test_rade-4p ${test_launch-4p})
set_tests_properties(test_rade-4p PROPERTIES PROCESSORS 4)
set_property(TEST test_rade-4p PROPERTY LABELS UNIT RADIATION)
