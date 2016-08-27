cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake"
                      "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules"
                      "${CMAKE_CURRENT_LIST_DIR}/../cmake/Utils")

include(NoInSourceBuild)
include(ExternalProject)

project(TruchasTPL C CXX Fortran)

# Fix a few things that cmake gets wrong for the NAG compiler.
if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  set(CMAKE_Fortran_COMPILE_OPTIONS_PIE "-PIC")
  include(NAGFortranCompilerVersion)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install"
      CACHE PATH "TPL install prefix" FORCE)
  list(APPEND CMAKE_SYSTEM_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}")
endif()

option(BUILD_SHARED_LIBS "Find/build shared TPLs" ON)
option(ENABLE_MPI "Find/build MPI-parallel TPLs where appropriate" ON)
option(SEARCH_FOR_HDF5 "Search for suitable HDF5 library" ON)
option(SEARCH_FOR_NETCDF "Search for suitable NetCDF library" ON)
option(SEARCH_FOR_EXODUS "Search for suitable Exodus library" ON)
option(SEARCH_FOR_HYPRE "Search for suitable HYPRE library" ON)
option(SEARCH_FOR_HYPRE "Search for suitable HYPRE library" ON)
option(SEARCH_FOR_YAJL "Search for suitable YAJL library" ON)
option(SEARCH_FOR_PETACA "Search for suitable Petaca library" ON)
option(SEARCH_FOR_SWIG "Search for suitable swig executable" ON)

# This really doesn't do anything now.  Originally, only .so were
# searched for when shared, and only .a when not shared.  This does not
# work for static (not shared) because some low-level system libraries
# like libm are only available as .so.  This would trip up the FindHDF
# module badly.  It boils down to what we mean by "shared" and "static"
# and how to enforce it.  Shared probably means, make shared libraries
# but use existing static ones if found, and static either the converse,
# or something more restrictive like build static libraries but use
# existing shared libraries if they are not one of the ones we can build?
if(BUILD_SHARED_LIBS)
  #message(STATUS "Finding and building shared libraries only.")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX}
                                  ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
  #message(STATUS "Finding and building static libraries only.")
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX}
                                  ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

# This macro converts a cmake-style list of library paths (full) to a string
# of -L and and corresponding -l arguments suitable for a link line.
macro(lib_list_to_link_arg var)
  set(${var})
  foreach(lib ${ARGN})
    get_filename_component(libdir ${lib} DIRECTORY)
    list(APPEND ${var} "-L${libdir}")
    get_filename_component(lib ${lib} NAME_WE)
    string(REGEX REPLACE "^lib" "" lib ${lib})
    list(APPEND ${var} "-l${lib}")
  endforeach()
  string(REPLACE ";" " " ${var} "${${var}}")
  unset(lib)
  unset(libdir)
endmacro()

# This macro converts a cmake list of directory paths to a string of -I arguments.
macro(dir_list_to_includes var)
  string(REGEX REPLACE "([^;]+)" "-I\\1" ${var} "${ARGN}")
  string(REPLACE ";" " " ${var} "${${var}}")
endmacro()

# This macro converts a cmake list of library directory paths to a string of -L arguments.
macro(dir_list_to_link_dirs var)
  string(REGEX REPLACE "([^;]+)" "-L\\1" ${var} "${ARGN}")
  string(REPLACE ";" " " ${var} "${${var}}")
endmacro()

# MPI (needed by HYPRE, HDF5)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# Python (needed by SWIG)
find_package(PythonInterp REQUIRED)

# Gather the C and Fortran flags that would be used by cmake into a single
# variable.  These will be passed to external project configure scripts.
set(cflags "${CMAKE_C_FLAGS}")
set(fflags "${CMAKE_Fortran_FLAGS}")
set(cppflags "${CMAKE_CXX_FLAGS}")
if(CMAKE_BUILD_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type)
  set(cflags "${cflags} ${CMAKE_C_FLAGS_${build_type}}")
  set(fflags "${fflags} ${CMAKE_Fortran_FLAGS_${build_type}}")
  set(cppflags "${cppflags} ${CMAKE_CXX_FLAGS_${build_type}}")
  string(STRIP "${cflags}" cflags)
  string(STRIP "${fflags}" fflags)
  string(STRIP "${cppflags}" cppflags)
endif()

# Most recent gcc compilers complain about HDF5 code
include(CheckCCompilerFlag)
check_c_compiler_flag("-std=c99" has_flag)
if(has_flag)
  set(cflags "${cflags} -std=c99")
endif()

#### Configure the TPLs as external projects as needed #########################

# Add the installation prefix to the list of paths that the various
# find_* commands will use when searching.
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}")

set(TARFILE_DIR ${PROJECT_SOURCE_DIR}/tarfiles)

set(projects_found)
set(projects_to_build)

include(hdf5_project)
include(netcdf_project)
include(exodus_project)
include(hypre_project)
include(yajl_project)
include(petaca_project)
include(swig_project)

# Print a summary of the configuration.
message(STATUS "")
string(REPLACE ";" " " project_list "${projects_found}")
if(projects_found)
  message(STATUS "Packages found: ${project_list}")
else()
  message(STATUS "No packages were found.")
endif()
string(REPLACE ";" " " project_list "${projects_to_build}")
if(projects_to_build)
  message(STATUS "Packages to build: ${project_list}")
  message(STATUS "")
  message(STATUS "Packages will be installed in ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "No packages to build; all were found.")
endif()
message(STATUS "")
