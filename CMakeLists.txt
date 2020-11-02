cmake_minimum_required(VERSION 3.16)

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake"
                      "${CMAKE_SOURCE_DIR}/cmake/Modules"
                      "${CMAKE_SOURCE_DIR}/cmake/Utils")

include(NoInSourceBuild)

project(Truchas LANGUAGES C CXX Fortran)

# Fix a few things that cmake gets wrong for the NAG compiler.
if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  set(CMAKE_Fortran_COMPILE_OPTIONS_PIE "-PIC")
  include(NAGFortranCompilerVersion)
endif()

# Create all the Fortran module files in a directory at the root build directory
set(Truchas_MODULE_DIR ${Truchas_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${Truchas_MODULE_DIR})

# ---------------------------------------------------------------------------- #
# Build Options
# ---------------------------------------------------------------------------- #

# Validate a user-supplied value for CMAKE_BUILD_TYPE.
if(CMAKE_BUILD_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)
  if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "^(DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be {Debug|Release|RelWithDebInfo|MinSizeRel}")
  endif()
endif()

# Build with support for dynamic loading of user functions
option(ENABLE_DYNAMIC_LOADING "Build with support for dynamic loading of user functions" True)

# Build with PIC (position in code)
option(WITH_PIC "Build with PIC support" True)
if(WITH_PIC)
  set(CMAKE_POSITION_INDEPENDENT_CODE True)
endif()

# Build shared libs
option(ENABLE_SHARED "Build shared libraries" True)
if(ENABLE_SHARED)
  set(BUILD_SHARED_LIBS True)
else()
  set(BUILD_SHARED_LIBS False)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

# Strip RPATH from installed binaries/libraries
option(ENABLE_INSTALL_RPATH "Use RPATH in installed binaries and libraries" True)

# ---------------------------------------------------------------------------- #
# External Dependencies
# ---------------------------------------------------------------------------- #

# MPI
find_package(MPI REQUIRED)

# Python and packages
include(FindPythonPackage)
find_package(PythonInterp 3.5 REQUIRED)
find_python_package(numpy 1.12)
find_python_package(h5py 2.6.0)

# BLAS and LAPACK -- these are not required
find_package(LAPACK)
find_package(BLAS)

# Rest are external dependencies that can be provided by the Truchas TPL build

set(CMAKE_PREFIX_PATH ${TRUCHAS_TPL_DIR})

# HYPRE library
find_package(HYPRE 2.15.1)
if(HYPRE_FOUND)
  if(NOT HYPRE_IS_PARALLEL)
    set(HYPRE_FOUND False)
    message(STATUS "Require parallel HYPRE library but found unsuitable serial library")
  endif()
endif()
if(NOT HYPRE_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Hypre installation")
endif()


# HDF5 library
set(HDF5_PREFER_PARALLEL True)
find_package(HDF5 "1.10.6" COMPONENTS C HL)
if(HDF5_FOUND)
  if(NOT HDF5_IS_PARALLEL)
    set(HDF5_FOUND False)
    message(STATUS "Require parallel HDF5 library but found unsuitable serial library")
  endif()
endif()
if(HDF5_FOUND)
  add_library(hdf5 INTERFACE)
  if(HDF5_C_LIBRARY)
    # hdf5 build from cmake
    target_link_libraries(hdf5 INTERFACE ${HDF5_C_HL_LIBRARY} ${HDF5_C_LIBRARY})
  else()
    # hdf5 build from autotools
    target_link_libraries(hdf5 INTERFACE ${HDF5_LIBRARIES})
  endif()
  target_include_directories(hdf5 INTERFACE ${HDF5_INCLUDE_DIRS})
  target_compile_definitions(hdf5 INTERFACE ${HDF5_DEFINITIONS})
else()
  message(FATAL_ERROR "Could NOT find a suitable HDF5 installation")
endif()

# NetCDF library
find_package(NetCDF "4.1.3")
if(NETCDF_FOUND)
  if(NOT NETCDF_HAS_NC4)
    message(STATUS "Found unsuitable NetCDF without the required netcdf-4 feature")
    set(NETCDF_FOUND False)
  endif()
endif()
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable NetCDF installation")
endif()

# ExodusII library
find_package(Exodus "514")
if(NOT EXODUS_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Exodus installation")
endif()

# Petaca library and its dependency libyajl
find_package(PETACA)
if(NOT PETACA_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Petaca installation")
endif()

# Chaco library
find_package(Chaco)
if(NOT Chaco_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Chaco installation")
endif()

# Chaparral library
find_package(Chaparral)
if(NOT Chaparral_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Chaparral installation")
endif()

# Scorpio library
find_package(Scorpio)
if(NOT Scorpio_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Scorpio installation")
endif()

# Portage library
option(USE_PORTAGE "Build with support for the Portage data mapping tool" OFF)
if (USE_PORTAGE)
  find_package(Portage)
  if(NOT PORTAGE_FOUND)
    message(FATAL_ERROR "Could NOT find a suitable Portage installation")
  endif()
endif()

# optional for building version information
find_package(Git)

# Test - need a swith here to enable testing
enable_testing()

# The value of SERIAL_COMMAND_PREFACE will preface the command line of serial
# unit tests.  Normally this is empty, but in the Cray Linux environment.
if(CRAY)
  set(SERIAL_COMMAND_PREFACE ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1)
else()
  set(SERIAL_COMMAND_PREFACE)
endif()

# ---------------------------------------------------------------------------- #
# Preprocessor definitions
# ----------------------------------------------------------------------------

# Add the compiler type to the global defines and the preprocessor
string(TOUPPER "${CMAKE_Fortran_COMPILER_ID}_COMPILER" Truchas_Compiler_Define)
add_compile_definitions(${Truchas_Compiler_Define} ${Truchas_Compiler_Define}_WORKAROUND)

# Specific compiler bug workarounds and defines. These are used by
# the custom Fortran preprocessing step (that need to be eliminated).
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "18.0")
    add_compile_definitions(NO_2008_FINDLOC)
  endif()
  add_compile_definitions(NO_2008_LHS_POLY_REALLOC)
  add_compile_definitions(INTEL_DPD200237439)
  add_compile_definitions(INTEL_BUG20180222)
  add_compile_definitions(INTEL_BUG20191228)
  add_compile_definitions(INTEL_BUG20191229)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "18.0.0")
    # This applies to 18.0.0 because of the version tweak number
    add_compile_definitions(INTEL18_WORKAROUND)
  endif()
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "19.0.0")
    add_compile_definitions(INTEL_SRN04330341)
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "7.0")
    add_compile_definitions(NO_2008_FINDLOC)
    add_compile_definitions(NO_2008_C_LOC)
  endif()
endif()

if(ENABLE_DYNAMIC_LOADING)
  add_compile_definitions(ENABLE_DYNAMIC_LOADING)
endif()

if(USE_PORTAGE)
  add_compile_definitions(USE_PORTAGE)
endif()

# More hardcoded defines - only used in drivers.f90 need to move
string(TIMESTAMP Truchas_BUILD_DATE "%Y-%m-%d %H:%M:%S")
site_name(Truchas_HOSTNAME)
# Flags for extra defines in drivers/drivers.F90
set(Truchas_INFO_FLAGS
  "BUILD_DATE='${Truchas_BUILD_DATE}'"
  "HOST_NAME='${Truchas_HOSTNAME}'"
  "ARCHITECTURE='${CMAKE_SYSTEM_PROCESSOR}'")

if(CMAKE_BUILD_TYPE)
  list(APPEND Truchas_INFO_FLAGS
    "COMPILER_FLAGS='${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}'")
else()
  list(APPEND Truchas_INFO_FLAGS "COMPILER_FLAGS='${CMAKE_Fortran_FLAGS}'")
endif()

# ---------------------------------------------------------------------------- #
# Installation Definitions
# ---------------------------------------------------------------------------- #

# Need the paths defined before adding src to the build

# Define the build id
include(GetTruchasBuildId)
get_truchas_build_id(Truchas_BUILD_ID)

# Default installation location is the Truchas_SOURCE_DIR
set(Truchas_DFLT_INSTALL_PREFIX ${Truchas_SOURCE_DIR}/install/${Truchas_BUILD_ID})
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${Truchas_DFLT_INSTALL_PREFIX}
      CACHE PATH "Installation prefix" FORCE)
endif()

# Binary installation path, override with -D Truchas_BIN_INSTALL_DIR
if (NOT Truchas_BIN_INSTALL_DIR)
  set(Truchas_BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
endif()

# Library installation path, override with -D Truchas_LIBRARY_INSTALL_DIR
if (NOT Truchas_LIBRARY_INSTALL_DIR)
  set(Truchas_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Module installation path, override with -D Truchas_MODULE_INSTALL_DIR
if (NOT Truchas_MODULE_INSTALL_DIR)
  set(Truchas_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Copy entire include directory to  Truchas_MODULE_INSTALL_DIR
install(CODE "MESSAGE(STATUS \"Installing: modules in ${Truchas_MODULE_INSTALL_DIR}\")")
install(CODE
        "EXECUTE_PROCESS ( COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${Truchas_MODULE_DIR}\" \"${Truchas_MODULE_INSTALL_DIR}\")"
)

# Configuration (CMake) file installation location
if (NOT Truchas_CONFIG_INSTALL_DIR)
  set(Truchas_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/truchas)
endif()

# ---------------------------------------------------------------------------- #
# RPATH handling  (http://www.cmake.org/Wiki/CMake_RPATH_handling)
# ---------------------------------------------------------------------------- #

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Library install directories
set(Truchas_LIBRARY_INSTALL_DIRS
     ${Truchas_LIBRARY_INSTALL_DIR}
     ${TRUCHAS_TPL_DIR}/lib)

# Add RPATH to binaries,libraries, *so, etc.
if (ENABLE_INSTALL_RPATH)

  # Only want non-system paths
  foreach(dir ${Truchas_LIBRARY_INSTALL_DIRS})
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${dir}" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
       list(APPEND CMAKE_INSTALL_RPATH "${dir}")
    endif("${isSystemDir}" STREQUAL "-1")
  endforeach()

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

endif()

# ---------------------------------------------------------------------------- #
# Source Directories
# ---------------------------------------------------------------------------- #

add_subdirectory(src)
add_subdirectory(test)

option(BUILD_HTML "Build HTML documentation (requires Sphinx)" OFF)
if(BUILD_HTML)
  find_package(Sphinx REQUIRED)
  add_subdirectory(doc/sphinx)
endif()

# ---------------------------------------------------------------------------- #
# Export Definitions
# ---------------------------------------------------------------------------- #

# Create the target file
#set(Truchas_TARGETS_FILE truchas-targets.cmake)
#export(truchas truchas_exe FILE ${Truchas_TARGETS_FILE})

# Create truchas-config.cmake (find_package(Truchas))
#set(Truchas_CONFIG_FILE ${Truchas_BINARY_DIR}/truchas-config.cmake)
#configure_file(${Truchas_SOURCE_DIR}/truchas-config.in
#               ${Truchas_CONFIG_FILE}
#	       @ONLY)
#install(FILES ${Truchas_CONFIG_FILE} ${Truchas_TARGETS_FILE}
#        DESTINATION ${Truchas_CONFIG_INSTALL_DIR})

# Print a summary of the build configuration
include(configuration_summary)
