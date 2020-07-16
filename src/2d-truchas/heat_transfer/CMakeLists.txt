#TODO: finalize build system
project(HeatTransfer2D Fortran)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# --- Library ht2d
set(HT2D_SOURCE_FILES
    mfd_2d_disc_type.F90
    mfd_2d_diff_matrix_type.F90
    mfd_2d_diff_precon_type.F90
    HT_2d_model_type.F90
    HT_2d_norm_type.F90
    HT_2d_precon_type.F90
    HT_2d_idaesol_model_type.F90
    HT_2d_solver_type.F90)

add_library(ht2d ${HT2D_SOURCE_FILES})
target_link_libraries(ht2d truchas ${Truchas_LINK_LIBS})

# --- ctest binaries
add_subdirectory(Test)

