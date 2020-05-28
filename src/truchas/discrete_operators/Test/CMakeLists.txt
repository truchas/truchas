project(DiscOpTest Fortran)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

add_executable(test_cell_grad_type test_cell_grad_type.F90)
target_link_libraries(test_cell_grad_type truchas ${Truchas_LINK_LIBS})

# 1-processor cell_grad_type test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./test_cell_grad_type ${DiscOpTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(cell_grad_type-1p ${test_launch})
set_tests_properties(cell_grad_type-1p PROPERTIES PROCESSORS 1)
set_property(TEST cell_grad_type-1p PROPERTY LABELS UNIT)

# 4-processor cell_grad_type test
set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test_cell_grad_type ${DiscOpTest_SOURCE_DIR} ${MPIEXEC_POSTFLAGS})
add_test(cell_grad_type-4p ${test_launch})
set_tests_properties(cell_grad_type-4p PROPERTIES PROCESSORS 4)
set_property(TEST cell_grad_type-4p PROPERTY LABELS UNIT)

# test development incomplete -- not working
#add_executable(test_linear_consistency test_linear_consistency.F90)
#target_link_libraries(test_linear_consistency truchas ${Truchas_LINK_LIBS})
