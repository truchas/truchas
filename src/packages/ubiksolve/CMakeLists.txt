# ############################################################################ #
#                                                                              #
# UbikSolve (CMake)                                                               #
#                                                                              #
# ############################################################################ #


# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(UbikSolve Fortran)

# Version
set(UbikSolve_MAJOR_VERSION 2)
set(UbikSolve_MINOR_VERSION 3)
set(UbikSolve_PATCH_VERSION b)
set(UbikSolve_VERSION ${UbikSolve_MAJOR_VERSION}.${UbikSolve_MINOR_VERSION}.${UbikSolve_PATCH_VERSION})

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${UbikSolve_SOURCE_DIR}")

# Create all the Fortran module files in a directory at the root build directory
set(UbikSolve_MODULE_DIR ${UbikSolve_BINARY_DIR}/modules)
global_set(UbikSolve_MODULE_DIR ${UbikSolve_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${UbikSolve_MODULE_DIR})


# ---------------------------------------------------------------------------- #
# Build Options 
# ---------------------------------------------------------------------------- #
include(CMakeDependentOption)

# PGSLib
set(ENABLE_PGSLib True CACHE BOOL "Build UbikSolve with PGSLib" FORCE)
if ( ENABLE_MPI )
  set(PGSLib pgslib_mpi)
else()
  set(PGSLib pgslib_serial)
endif()

# Tests
cmake_dependent_option(ENABLE_TESTS "Build the test suite" False
                       "NOT ENABLE_PGSLib" False)
if (ENABLE_TESTS)
  enable_testing()
endif()

# ---------------------------------------------------------------------------- #
# Preprocessor definitions 
# ---------------------------------------------------------------------------- #

# Preprocessor - allow this to be a list to define global flags easily
if ( NOT UbikSolve_PREPROCESSOR )
  set(UbikSolve_PREPROCESSOR cpp -ansi)
endif()

# Add the compiler type to the global defines and the preprocessor
string(TOUPPER "-D${CMAKE_Fortran_COMPILER_ID}_COMPILER" UbikSolve_Compiler_Define)
add_definitions(${UbikSolve_Compiler_Define})
list(APPEND UbikSolve_PREPROCESSOR ${UbikSolve_Compiler_Define})

# Add any user defined preprocessor flags
if (UbikSolve_FPP_FLAGS)
  list(APPEND UbikSolve_PREPROCESSOR ${UbikSolve_FPP_FLAGS})
endif()

# Add the PGSLIb defines
if (ENABLE_PGSLib)
  list(APPEND UbikSolve_PREPROCESSOR -Duse_PGSLib)
endif()

# ---------------------------------------------------------------------------- #
# Installation Definitions 
# ---------------------------------------------------------------------------- #

# Need the paths defined before adding src to the build

# Library installation path, override with -D UbikSolve_LIBRARY_INSTALL_DIR
if (NOT UbikSolve_LIBRARY_INSTALL_DIR)
  set(UbikSolve_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Module installation path, override with -D UbikSolve_MODULE_INSTALL_DIR
if (NOT UbikSolve_MODULE_INSTALL_DIR)
  set(UbikSolve_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Copy entire include directory to  UbikSolve_MODULE_INSTALL_DIR
install(CODE "MESSAGE(STATUS \"Installing: modules in ${UbikSolve_MODULE_INSTALL_DIR}\")")
install(CODE
        "EXECUTE_PROCESS ( COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${UbikSolve_MODULE_DIR}\" \"${UbikSolve_MODULE_INSTALL_DIR}\")"
)

if (NOT UbikSolve_CONFIG_INSTALL_DIR)
  set(UbikSolve_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/ubiksolve)
endif()

# ---------------------------------------------------------------------------- #
# Source Directories 
# ---------------------------------------------------------------------------- #

add_subdirectory(src)

if (ENABLE_TESTS)
  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------- #
# Export Definitions 
# ---------------------------------------------------------------------------- #

# Create the target file
set(UbikSolve_TARGETS_FILE ubiksolve-targets.cmake)
install(EXPORT ubiksolve
        FILE ${UbikSolve_TARGETS_FILE}
        DESTINATION ${UbikSolve_CONFIG_INSTALL_DIR})

# Create ubiksolve-config.cmake (find_package(UbikSolve))
set(UbikSolve_CONFIG_FILE ${UbikSolve_BINARY_DIR}/ubiksolve-config.cmake)
configure_file(${UbikSolve_SOURCE_DIR}/ubiksolve-config.in
               ${UbikSolve_CONFIG_FILE}
	       @ONLY)
install(FILES ${UbikSolve_CONFIG_FILE} 
        DESTINATION ${UbikSolve_CONFIG_INSTALL_DIR})


