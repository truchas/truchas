project(TruchasToolsPP)

message(STATUS "configuring the Truchas post_processors directory")

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# Add the Danu module directory
include_directories(${Danu_Fortran_MODULE_DIR})

# --- Preprocessing flags
set(fpp_flags_list 
    ${Truchas_FPP_FLAGS} -I${TruchasToolsPP_SOURCE_DIR})
list(APPEND fpp_flags_list -I${Danu_Fortran_MODULE_DIR})
build_whitespace_string(fpp_flags ${fpp_flags_list})

# --- write_restart binary
fortran_preprocess_files(WRITE_RESTART_FILES
                         FILES
			      f90_assert.F90
			      h5out_type.F90
			      mapped_restart.F90
			      write_restart.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessWriteRestartFiles)
set(fc_flags -I${NETCDF_INCLUDE_DIR})
list(APPEND fc_flags -I${Danu_Fortran_MODULE_DIR})
if(Fortran_COMPILER_IS_INTEL)
  list(APPEND fc_flags "-assume realloc_lhs")
endif()
build_whitespace_string(WRITE_RESTART_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${WRITE_RESTART_FILES} PROPERTIES
                            COMPILE_FLAGS ${WRITE_RESTART_COMPILE_FLAGS})
add_executable(write_restart ${WRITE_RESTART_FILES})
add_dependencies(write_restart ProcessWriteRestartFiles)
target_link_libraries(write_restart truchas ${Truchas_LINK_LIBS})

# --- write_probes binary
fortran_preprocess_files(WRITE_PROBES_FILES
                         FILES
			      f90_assert.F90
			      h5out_type.F90
			      write_probes.F90
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${fpp_flags_list}
			 PROCESS_TARGET ProcessWriteProbeFiles)
set_source_files_properties(${WRITE_PROBE_FILES} PROPERTIES
                            COMPILE_FLAGS ${WRITE_RESTART_COMPILE_FLAGS})
add_executable(write_probe ${WRITE_PROBES_FILES})
add_dependencies(write_probe ProcessWriteProbeFiles)
target_link_libraries(write_probe truchas ${Truchas_LINK_LIBS})


install(TARGETS write_restart write_probe
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})

