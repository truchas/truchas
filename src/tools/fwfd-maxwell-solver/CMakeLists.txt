project(FWFD-Maxwell-Solver Fortran)

# Want modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

# --- driver program binary
add_executable(fwfd_solver fdme_sim_type.F90 main.F90)
target_link_libraries(fwfd_solver truchas ${Truchas_LINK_LIBS})
install(TARGETS fwfd_solver
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
