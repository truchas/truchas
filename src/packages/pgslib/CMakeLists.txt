# ############################################################################ #
#                                                                              #
# PGSLib (CMake)                                                               #
#                                                                              #
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(PGSLib C Fortran)

# Version
set(PGSLib_MAJOR_VERSION 3)
set(PGSLib_MINOR_VERSION 3)
set(PGSLib_PATCH_VERSION 0)
set(PGSLib_VERSION ${PGSLib_MAJOR_VERSION}.${PGSLib_MINOR_VERSION}.${PGSLib_PATCH_VERSION})

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PGSLib_SOURCE_DIR}")

# Create all the Fortran module files in a directory at the root build directory
set(PGSLib_MODULE_DIR ${PGSLib_BINARY_DIR}/include CACHE INTERNAL "")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PGSLib_MODULE_DIR})

# ---------------------------------------------------------------------------- #
# Preprocessor definitions
# ---------------------------------------------------------------------------- #

# Add the compiler type to the global defines and the preprocessor
string(TOUPPER "-D${CMAKE_Fortran_COMPILER_ID}_COMPILER" PGSLib_Compiler_Define)
add_definitions(${PGSLib_Compiler_Define})

# ---------------------------------------------------------------------------- #
# Installation Definitions
# ---------------------------------------------------------------------------- #

# Need the paths defined before adding src to the build

# Library installation path, override with -D PGSLib_LIBRARY_INSTALL_DIR
if (NOT PGSLib_LIBRARY_INSTALL_DIR)
  set(PGSLib_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Module installation path, override with -D PGSLib_MODULE_INSTALL_DIR
if (NOT PGSLib_MODULE_INSTALL_DIR)
  set(PGSLib_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Copy entire include directory to  PGSLib_MODULE_INSTALL_DIR
install(CODE "MESSAGE(STATUS \"Installing: modules in ${PGSLib_MODULE_INSTALL_DIR}\")")
install(CODE
        "EXECUTE_PROCESS ( COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${PGSLib_MODULE_DIR}\" \"${PGSLib_MODULE_INSTALL_DIR}\")"
)

if (NOT PGSLib_CONFIG_INSTALL_DIR)
  set(PGSLib_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/pgslib)
endif()

# ---------------------------------------------------------------------------- #
# Source Directories
# ---------------------------------------------------------------------------- #

add_subdirectory(src)

# ---------------------------------------------------------------------------- #
# Export Definitions
# ---------------------------------------------------------------------------- #

# Create the target file
set(PGSLib_TARGETS_FILE pgslib-targets.cmake)
set(export_targets pgslib_fe)
if (ENABLE_MPI)
  list(APPEND export_targets pgslib_mpi)
endif()
foreach( t ${export_targets} )
  install(EXPORT ${t}
          DESTINATION ${PGSLib_CONFIG_INSTALL_DIR}
          FILE ${PGSLib_TARGETS_FILE})
endforeach()


# Create pgslib-config.cmake (find_package(PGSLib))
set(PGSLib_CONFIG_FILE ${PGSLib_BINARY_DIR}/pgslib-config.cmake)
configure_file(${PGSLib_SOURCE_DIR}/pgslib-config.in
               ${PGSLib_CONFIG_FILE}
	       @ONLY)
install(FILES ${PGSLib_CONFIG_FILE}
        DESTINATION ${PGSLib_CONFIG_INSTALL_DIR})
