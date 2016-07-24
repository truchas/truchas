project(TruchasToolsScramble Fortran)

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# --- Preprocessor
# To avoid overwriting VERSION
set(SCRAMBLE_PREPROCESSOR cpp -traditional-cpp)


# --- Scramble
fortran_preprocess_files(SCRAMBLE_SOURCE_FILES
                         FILES 
			     f90_assert.F90
			     scramble_proc.F90 
			     command_line.F90
			     scramble.F90
			 FPP_EXECUTABLE ${SCRAMBLE_PREPROCESSOR}
			 FPP_FLAGS -I${TruchasToolsScramble_SOURCE_DIR}
			 PROCESS_TARGET ProcessScrambleFiles)

add_executable(scramble ${SCRAMBLE_SOURCE_FILES})
add_dependencies(scramble ProcessScrambleFiles)
target_link_libraries(scramble truchas ${Truchas_LINK_LIBS})

# --- Install rules 		    
install(TARGETS scramble
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})

