project(DiscOpTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${Danu_Fortran_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

add_executable(test_cell_grad_type test_cell_grad_type.F90)
target_link_libraries(test_cell_grad_type truchas ${Truchas_LINK_LIBS})

if(ENABLE_MPI)
  set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_cell_grad_type ${DiscOpTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
else()
  set(test_launch test_cell_grad_type ${DiscOpTest_SOURCE_DIR})
endif()
add_test(cell_grad_type ${test_launch})

if(ENABLE_MPI)
  set_tests_properties(cell_grad_type PROPERTIES PROCESSORS 4)
endif()
set_property(TEST cell_grad_type PROPERTY LABELS UNIT)

# test development incomplete -- not working
#add_executable(test_linear_consistency test_linear_consistency.F90)
#target_link_libraries(test_linear_consistency truchas ${Truchas_LINK_LIBS})
