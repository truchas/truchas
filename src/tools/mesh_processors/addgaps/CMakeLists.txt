project(TruchasToolsAddgaps Fortran)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# --- AddGaps
set(ADDGAPS_SOURCE_FILES
  f90_assert.F90
  addgaps_proc.F90
  command_line.F90
  addgaps.F90)

add_executable(addgaps ${ADDGAPS_SOURCE_FILES})
target_link_libraries(addgaps truchas)

# --- Install rules
install(TARGETS addgaps
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
