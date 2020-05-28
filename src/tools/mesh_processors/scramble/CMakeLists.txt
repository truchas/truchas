project(TruchasToolsScramble Fortran)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# --- Scramble
set(SCRAMBLE_SOURCE_FILES
  f90_assert.F90
  scramble_proc.F90
  command_line.F90
  scramble.F90)

add_executable(scramble ${SCRAMBLE_SOURCE_FILES})
target_link_libraries(scramble truchas ${Truchas_LINK_LIBS})

# --- Install rules
install(TARGETS scramble
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
