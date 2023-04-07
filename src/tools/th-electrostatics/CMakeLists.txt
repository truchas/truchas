project(TH-Electrostatics Fortran)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

# --- driver program binary
add_executable(prog thes_sim_type.F90 thes_solver_type.F90 thes_bc_type.F90 main.F90)
target_link_libraries(prog truchas ${Truchas_LINK_LIBS})
set_target_properties(prog PROPERTIES OUTPUT_NAME thes_solver)
install(TARGETS prog
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
