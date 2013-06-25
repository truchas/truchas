project(TruchasToolsMesh)

message(STATUS "configuring the Truchas post_processors directory")

# --- CMake modules
include(FortranPreprocessFiles)

# --- Include paths

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# Truchas module directory
include_directories(${Truchas_MODULE_DIR})

# Add the Danu module directory
if(ENABLE_Danu)
  include_directories(${Danu_Fortran_MODULE_DIR})
endif()

# --- Preprocessing flags
set(fpp_flags_list 
    ${Truchas_FPP_FLAGS} -I${TruchasToolsMesh_SOURCE_DIR})
if(ENABLE_Danu)
  list(APPEND fpp_flags_list -I${Danu_Fortran_MODULE_DIR})
endif()
build_whitespace_string(fpp_flags ${fpp_flags_list})


# --- mesh_translate binary
add_executable(meshtranslate meshtranslate.f90)


# --- Add subdirectories 

# addgaps 
add_subdirectory(addgaps)

# --- Install rules 		    
install(TARGETS meshtranslate
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
	ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})


