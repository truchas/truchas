# ############################################################################ #
#                                                                              #
# Chaco  (CMake)                                                               #
#                                                                              #
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(Chaco C)

# Version
set(Chaco_MAJOR_VERSION 2)
set(Chaco_MINOR_VERSION 2)
set(Chaco_PATCH_VERSION 1)
set(Chaco_VERSION ${Chaco_MAJOR_VERSION}.${Chaco_MINOR_VERSION}.${Chaco_PATCH_VERSION})

# ---------------------------------------------------------------------------- #
# Installation definitions 
# ---------------------------------------------------------------------------- #

# Library files
if ( NOT Chaco_LIBRARY_INSTALL_DIR )
  set(Chaco_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# CMake Config files
if ( NOT Chaco_CONFIG_INSTALL_DIR )
  set(Chaco_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/chaco)
endif()  

# ---------------------------------------------------------------------------- #
# Main source directory
# ---------------------------------------------------------------------------- #

add_subdirectory(src)

# ---------------------------------------------------------------------------- #
# Export Definitions 
# ---------------------------------------------------------------------------- #

# Create the target file
set(Chaco_TARGETS_FILE chaco-targets.cmake)
install(EXPORT chaco 
        DESTINATION ${Chaco_CONFIG_INSTALL_DIR}
	FILE ${Chaco_TARGETS_FILE})

# Create chaco-config.cmake (find_package(Chaco))
set(Chaco_CONFIG_FILE ${Chaco_BINARY_DIR}/chaco-config.cmake)
configure_file(${Chaco_SOURCE_DIR}/chaco-config.in
               ${Chaco_CONFIG_FILE}
	       @ONLY)
install(FILES ${Chaco_CONFIG_FILE}
        DESTINATION ${Chaco_CONFIG_INSTALL_DIR})
