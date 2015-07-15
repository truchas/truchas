project(MFDdiscTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})

add_executable(test_mfd_polygon_matrix_inv test_mfd_polygon_matrix_inv.F90)
target_link_libraries(test_mfd_polygon_matrix_inv truchas ${Truchas_LINK_LIBS})

set(test_launch test_mfd_polygon_matrix_inv)

add_test(mfd_polygon_matrix_inv ${test_launch})

#set_tests_properties(hypre_cell_scheme_matrix_type PROPERTIES PROCESSORS 4)
#set_property(TEST hypre_cell_scheme_matrix_type PROPERTY LABELS UNIT)