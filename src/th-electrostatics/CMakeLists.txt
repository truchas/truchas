project(TH-Electrostatics Fortran)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

## --- Library ht2d
#set(HT2D_SOURCE_FILES
#    mfd_2d_disc_type.F90
#    mfd_2d_diff_matrix_type.F90
#    mfd_2d_diff_precon_type.F90
#    HT_2d_model_type.F90
#    HT_2d_norm_type.F90
#    HT_2d_precon_type.F90
#    HT_2d_idaesol_model_type.F90
#    HT_2d_solver_type.F90
#    HT_2d_sim_type.F90
#    new_TofH_type.F90
#)
#
#add_library(ht2d ${HT2D_SOURCE_FILES})
#target_link_libraries(ht2d truchas truchas-extra ${Truchas_LINK_LIBS})

# --- driver program binary
add_executable(prog th_electrostatics_sim_type.F90 main.F90)
target_link_libraries(prog truchas ${Truchas_LINK_LIBS})
set_target_properties(prog PROPERTIES OUTPUT_NAME th-electrostatics)
