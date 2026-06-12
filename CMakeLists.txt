cmake_minimum_required(VERSION 3.20.2)

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake"
                      "${CMAKE_SOURCE_DIR}/cmake/Modules"
                      "${CMAKE_SOURCE_DIR}/cmake/Utils")

include(NoInSourceBuild)

project(Truchas LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

option(ENABLE_STD_MOD_PROC_NAME "Build with -assume std_mod_proc_name when using Intel" OFF)

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

# LAPACK (and BLAS)
find_package(LAPACK REQUIRED)

# Rest are external dependencies that can be provided by the Truchas TPL build

set(CMAKE_PREFIX_PATH ${TRUCHAS_TPL_DIR})

# HYPRE library
find_package(HYPRE 2.15.1 CONFIG REQUIRED)
if(NOT HYPRE_WITH_MPI)
  message(FATAL_ERROR "Require parallel HYPRE library but found unsuitable serial library")
endif()

# HDF5 library
find_package(HDF5 CONFIG REQUIRED COMPONENTS C HL shared)
if(NOT HDF5_ENABLE_PARALLEL)
  message(FATAL_ERROR "Require parallel HDF5 library but found unsuitable serial library")
endif()
if(NOT TARGET hdf5::hdf5)
  add_library(hdf5::hdf5 INTERFACE IMPORTED)
  target_link_libraries(hdf5::hdf5 INTERFACE hdf5_hl-shared hdf5-shared)
endif()

# Use this form, and remove the fVTKHDF_MPI_SUPPORT check below, after fVTKHDF
# fixes component reporting in its next package config:
# find_package(fVTKHDF 0.6.0 CONFIG REQUIRED COMPONENTS MPI)
find_package(fVTKHDF 0.6.0 CONFIG REQUIRED)
if(NOT fVTKHDF_MPI_SUPPORT)
  message(FATAL_ERROR "Require fVTKHDF with MPI support")
endif()

# MUMPS library
option(USE_MUMPS "Build with support for the MUMPS direct solver library" OFF)
if (USE_MUMPS)
  find_package(MUMPS 5.7 REQUIRED COMPONENTS d z)
endif()

# NetCDF library
find_package(netCDF "4.1.3" CONFIG REQUIRED)
if(NOT netCDF_HAS_NC4)
  message(FATAL_ERROR "Require NetCDF with the netcdf-4 feature")
endif()

# ExodusII library
find_package(SEACASExodus CONFIG REQUIRED)

# Petaca library and its dependency libyajl
find_package(PETACA)
if(NOT PETACA_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Petaca installation")
endif()

# Chaparral library
find_package(Chaparral CONFIG REQUIRED)
if(NOT TARGET Chaparral::Chaparral AND TARGET vflib)
  add_library(Chaparral::Chaparral INTERFACE IMPORTED)
  target_link_libraries(Chaparral::Chaparral INTERFACE vflib)
endif()
if(NOT TARGET Chaparral::Chaparral)
  message(FATAL_ERROR "Chaparral package config did not define a usable library target")
endif()

# Scorpio library
find_package(Scorpio)
if(NOT Scorpio_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Scorpio installation")
endif()

# Portage library
option(USE_PORTAGE "Build with support for the Portage data mapping tool" OFF)
if (USE_PORTAGE)
  find_package(portage)
  if(NOT portage_FOUND)
    message(FATAL_ERROR "Could NOT find a suitable Portage installation")
  endif()
endif()

# Metis library
find_package(METIS 5.1.0)
if(NOT METIS_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Metis installation")
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


add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-standard-semantics>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-standard-semantics>"
  "$<$<AND:$<NOT:$<BOOL:${ENABLE_STD_MOD_PROC_NAME}>>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:SHELL:-assume nostd_mod_proc_name>"
  "$<$<AND:$<NOT:$<BOOL:${ENABLE_STD_MOD_PROC_NAME}>>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:SHELL:-assume nostd_mod_proc_name>"
)

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>)

# Add the compiler type to the global defines and the preprocessor
string(TOUPPER "${CMAKE_Fortran_COMPILER_ID}_COMPILER" Truchas_Compiler_Define)
add_compile_definitions(${Truchas_Compiler_Define} ${Truchas_Compiler_Define}_WORKAROUND)

# Compiler-specific bug workarounds and defines.
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel OR CMAKE_Fortran_COMPILER_ID MATCHES IntelLLVM)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "19.0.4")
    # Needed by any version earlier than 19.0.4
    add_compile_definitions(INTEL_BUG20180115)
  endif()
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "20.0")
    # Needed by any version earlier than oneAPI 2021.1
    add_compile_definitions(INTEL_BUG20191229)
    add_compile_definitions(INTEL_BUG20190903)
    add_compile_definitions(INTEL_BUG20200721)
    add_compile_definitions(INTEL_BUG20211122)
  endif()
  # Needed by any version through at least 2021.2
  add_compile_definitions(NO_2008_LHS_POLY_REALLOC)
  add_compile_definitions(INTEL_BUG20121013A)
  add_compile_definitions(INTEL_BUG20180222)
  add_compile_definitions(INTEL_BUG20191228)
  add_compile_definitions(INTEL_BUG20210619)
  add_compile_definitions(INTEL_BUG20240327)
  add_compile_definitions(INTEL_BUG20241231)
  add_compile_definitions(INTEL_BUG20250605)
  add_compile_definitions(INTEL_BUG20250709)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "7.1")
    add_compile_definitions(NAG_BUG20220103)
    add_compile_definitions(NAG_BUG20260521)
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  add_compile_definitions(GNU_PR49213)
  add_compile_definitions(GNU_PR69563)
  add_compile_definitions(GNU_PR86277)
  add_compile_definitions(GNU_PR93762)
  add_compile_definitions(GNU_PR98141)
  add_compile_definitions(GNU_ISSUE437)
  add_compile_definitions(GNU_PR103394)
  add_compile_definitions(GNU_PR117774)
  add_compile_definitions(GNU_PR119986)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.3.0")
    add_compile_definitions(NO_2008_LHS_POLY_REALLOC) # see issue 475
  endif()
endif()

if(APPLE)
  add_compile_definitions(AVOID_MPI_IN_PLACE)
endif()

if(ENABLE_DYNAMIC_LOADING)
  add_compile_definitions(ENABLE_DYNAMIC_LOADING)
endif()

if(USE_MUMPS)
  add_compile_definitions(USE_MUMPS)
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
  find_program(SPHINX_EXECUTABLE NAMES sphinx-build
    DOC "Path to sphinx-build executable")
  if(NOT SPHINX_EXECUTABLE)
    message(FATAL_ERROR "sphinx-build not found. Please install Sphinx.")
  endif()
  add_subdirectory(doc/sphinx)
  add_subdirectory(doc/reference-manual)
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
