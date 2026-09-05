#TODO: finalize build system
project(HeatTransfer2D Fortran)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# --- Library ht2d
set(HT2D_SOURCE_FILES
    # Temporary bridge until refactor-diffusion-solver provides these through
    # the truchas library.
    new_nka_type.F90
    new_state_history_type.F90
    new_idaesol_type.F90
    mfd_2d_disc_type.F90
    mfd_2d_diff_matrix_type.F90
    mfd_2d_diff_precon_type.F90
    ht_2d_vector_type.F90
    thermal_bc_factory_type.F90
    ht_2d_source_factory_type.F90
    ht_2d_model_type.F90
    ht_2d_norm_type.F90
    ht_2d_precon_type.F90
    ht_2d_idaesol_model_type.F90
    ht_2d_ic_solver_type.F90
    ht_2d_solver_type.F90
    ht_2d_sim_type.F90
    ht_2d_vtkhdf_writer_type.F90
    ht_2d_tofh_type.F90
)

add_library(ht2d ${HT2D_SOURCE_FILES})
set_target_properties(ht2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(ht2d PUBLIC truchas truchas-extra ${Truchas_LINK_LIBS})
target_include_directories(ht2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

# --- ctest binaries
add_subdirectory(test)
