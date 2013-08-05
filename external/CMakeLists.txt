# #  -*- mode: cmake -*-
# ############################################################################ #
#                                                                              #
# TruchasExternal                                                              #
#                                                                              #
# ############################################################################ #

project(TruchasExternal)

message(STATUS "configuring Truchas external projects")

# --- Add the local cmake directory 
list(APPEND CMAKE_MODULE_PATH ${TruchasExternal_SOURCE_DIR}/cmake)

# --- Define the external project installation  
if (NOT TruchasExternal_INSTALL_PREFIX)
  set(TruchasExternal_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/external
      CACHE PATH "Truchas external projects install prefix")
  set(TruchasExternal_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX} PARENT_SCOPE)  
endif()

# --- Useful directory variables
set(TruchasExternal_TEMPLATES_DIR ${TruchasExternal_SOURCE_DIR}/templates)
set(TruchasExternal_ARCHIVE_DIR   ${TruchasExternal_SOURCE_DIR}/tarfiles)
set(TruchasExternal_PATCHES_DIR   ${TruchasExternal_SOURCE_DIR}/patches)

# --- Search and verify the external packages
if ( ENABLE_ExternalSearch )
  include(Verify_HYPRE)
  include(Verify_SWIG)
  include(Verify_HDF5)
  include(Verify_NetCDF)
else(ENABLE_ExternalSearch)  
  set(HYPRE_VERIFIED False)
  set(SWIG_VERIFIED False)
  set(HDF5_VERIFIED False)
  set(NETCDF_VERIFIED False)
endif(ENABLE_ExternalSearch)  

# Set a flag to indicate that builds are required
include(BoolEval)
bool_eval(BUILD_EXTERNAL 
          NOT HYPRE_VERIFIED OR
	  NOT SWIG_VERIFIED OR
	  NOT HDF5_VERIFIED OR
	  NOT NETCDF_VERIFIED)

# --- External projects built for Truchas

if(BUILD_EXTERNAL)

  # Adds ExternalProject_Add to scope
  include(ExternalProject)

  # Create configure, build, install and test targets for each ExternalProject
  set_property(DIRECTORY PROPERTY
               EP_STEP_TARGETS download patch configure build install test)

  # CMake compiler settings for any package built with CMake
  set(TruchasExternal_CMAKE_COMPILER_ARGS
     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
     -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
     -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
     -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
     -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
     -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
     -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
     -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
     -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
     -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
     -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  )

  # CMake build settings
  set(TruchasExternal_CMAKE_BUILD_ARGS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  )

  # GNU AUTOCONF configure static/shared
  set(TruchasExternal_SHARED_SWITCH --disable-shared)
  if(BUILD_SHARED_LIBS)
    set(TruchasExternal_SHARED_SWITCH --enable-shared)
  endif()

  # Log all activity -- log files for each step in the timestamp directory
  set(TruchasExternal_LOG_OPTS
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_TEST      1
      LOG_INSTALL   1
  )

  # Loop through and add ExternalProject_* if not verified
  foreach(proj HDF5 NetCDF HYPRE SWIG)
    if ( NOT "${${proj}_VERIFIED}")
      include(ExternalProject_${proj})
    endif(NOT "${${proj}_VERIFIED}")
  endforeach()  

endif(BUILD_EXTERNAL)

# --- External project variables
# -   Other parts of the build will need these variables
# -   defined correctly. Varibles are either defined in Find*.cmake 
# -   or ExternalProject_*.cmake files. Without global settings
# -   these variables and their values are out of scope outside of
# -   this subdirectory.  

# ZLIB
global_set(ZLIB_FOUND          ${ZLIB_FOUND})
global_set(ZLIB_VERSION        ${ZLIB_VERSION})
global_set(ZLIB_VERSION_STRING ${ZLIB_VERSION_STRING})
global_set(ZLIB_LIBRARY        ${ZLIB_LIBRARY})
global_set(ZLIB_LIBRARIES      ${ZLIB_LIBRARIES})
global_set(ZLIB_INCLUDE_DIR    ${ZLIB_INCLUDE_DIR})
global_set(ZLIB_INCLUDE_DIRS   ${ZLIB_INCLUDE_DIRS})

# SWIG
global_set(SWIG_FOUND ${SWIG_FOUND})
global_set(SWIG_DIR ${SWIG_DIR})
global_set(SWIG_EXECUTABLE ${SWIG_EXECUTABLE})

# HYPRE
global_set(HYPRE_FOUND ${HYPRE_FOUND})
global_set(HYPRE_LIBRARY ${HYPRE_LIBRARY})
global_set(HYPRE_LIBRARIES ${HYPRE_LIBRARIES})
global_set(HYPRE_INCLUDE_DIR ${HYPRE_INCLUDE_DIR})
global_set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIRS})

# HDF5
global_set(HDF5_FOUND ${HDF5_FOUND})
global_set(HDF5_VERSION ${HDF5_VERSION})
global_set(HDF5_C_LIBRARY ${HDF5_C_LIBRARY})
global_set(HDF5_C_LIBRARIES ${HDF5_C_LIBRARIES})
global_set(HDF5_LINK_LIBRARIES ${HDF5_LINK_LIBRARIES})
global_set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR})
global_set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})

# NETCDF
global_set(NETCDF_FOUND ${NETCDF_FOUND})
global_set(NETCDF_C_LIBRARY ${NETCDF_C_LIBRARY})
global_set(NETCDF_C_LIBRARIES ${NETCDF_C_LIBRARIES})
global_set(NETCDF_Fortran_LIBRARY ${NETCDF_Fortran_LIBRARY})
global_set(NETCDF_Fortran_LIBRARIES ${NETCDF_Fortran_LIBRARIES})
global_set(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR})
global_set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIRS})
