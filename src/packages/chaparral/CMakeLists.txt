# ############################################################################ #
#                                                                              #
# Chaparral  (CMake)                                                               #
#                                                                              #
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(Chaparral C)

# Version
set(Chaparral_MAJOR_VERSION 3)
set(Chaparral_MINOR_VERSION 2)
set(Chaparral_VERSION ${Chaparral_MAJOR_VERSION}.${Chaparral_MINOR_VERSION})

# ---------------------------------------------------------------------------- #
# Build options
# ---------------------------------------------------------------------------- #

# Enable ExodusII
option(ENABLE_ExodusII "Build with Exodus" False)
if(ENABLE_ExodusII)
  message(FATAL_ERROR "ENABLE_ExodusII is not supported at this time")
endif()


# ---------------------------------------------------------------------------- #
# External dependencies
# ---------------------------------------------------------------------------- #

if (ENABLE_MPI)
  find_package(MPI COMPONENTS C)
endif()

if (ENABLE_ExodusII)
  find_package(ExodusII)
endif()

# ---------------------------------------------------------------------------- #
# Global defines
# ---------------------------------------------------------------------------- #
if(APPLE)
  add_definitions(-Dmacosx)
elseif(WIN32)
  add_definitions(-DWIN32)
elseif(UNIX)
  add_definitions(-Dlinux)
endif()

# ---------------------------------------------------------------------------- #
# Installation Definitions
# ---------------------------------------------------------------------------- #

# Need the paths defined before adding libraries and executables to the build

# Include installation path, override with -D Chaparral_INCLUDE_INSTALL_DIR
if (NOT Chaparral_INCLUDE_INSTALL_DIR)
  set(Chaparral_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Library installation path, override with -D Chaparral_LIBRARY_INSTALL_DIR
if (NOT Chaparral_LIBRARY_INSTALL_DIR)
  set(Chaparral_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Configuration installation path, override with -D Chaparral_CONFIG_INSTALL_DIR
if (NOT Chaparral_CONFIG_INSTALL_DIR)
  set(Chaparral_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/chaparral)
endif()

# ---------------------------------------------------------------------------- #
# Source directories
# ---------------------------------------------------------------------------- #

add_subdirectory(vflib)
add_subdirectory(tools)

# ---------------------------------------------------------------------------- #
# Export Definitions 
# ---------------------------------------------------------------------------- #

# Create the target file
set(Chaparral_TARGETS_FILE chaparral-targets.cmake)
install(EXPORT vflib
        DESTINATION ${Chaparral_CONFIG_INSTALL_DIR}
	FILE ${Chaparral_TARGETS_FILE})

# Create chaparral-config.cmake (find_package(Chaparral))
set(Chaparral_CONFIG_FILE ${Chaparral_BINARY_DIR}/chaparral-config.cmake)
configure_file(${Chaparral_SOURCE_DIR}/chaparral-config.in
               ${Chaparral_CONFIG_FILE}
	       @ONLY)
install(FILES ${Chaparral_CONFIG_FILE}
        DESTINATION ${Chaparral_CONFIG_INSTALL_DIR})


