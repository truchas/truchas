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
# -   or ExternalProject_*.cmake files. Without PARENT_SCOPE settings
# -   these variables and their values are out of scope outside of
# -   this subdirectory.  

# ZLIB
set(ZLIB_FOUND          ${ZLIB_FOUND} PARENT_SCOPE)
set(ZLIB_BUILD_TARGET   ${ZLIB_BUILD_TARGET} PARENT_SCOPE)
set(ZLIB_VERSION        ${ZLIB_VERSION} PARENT_SCOPE)
set(ZLIB_VERSION_STRING ${ZLIB_VERSION_STRING} PARENT_SCOPE)
set(ZLIB_LIBRARY        ${ZLIB_LIBRARY} PARENT_SCOPE)
set(ZLIB_LIBRARIES      ${ZLIB_LIBRARIES} PARENT_SCOPE)
set(ZLIB_INCLUDE_DIR    ${ZLIB_INCLUDE_DIR} PARENT_SCOPE)
set(ZLIB_INCLUDE_DIRS   ${ZLIB_INCLUDE_DIRS} PARENT_SCOPE)

# SWIG
set(SWIG_FOUND        ${SWIG_FOUND} PARENT_SCOPE)
set(SWIG_BUILD_TARGET ${SWIG_BUILD_TARGET} PARENT_SCOPE)
set(SWIG_DIR          ${SWIG_DIR} PARENT_SCOPE)
set(SWIG_EXECUTABLE   ${SWIG_EXECUTABLE} PARENT_SCOPE)

# HYPRE
set(HYPRE_FOUND        ${HYPRE_FOUND} PARENT_SCOPE)
set(HYPRE_BUILD_TARGET ${HYPRE_BUILD_TARGET} PARENT_SCOPE)
set(HYPRE_VERSION      ${HYPRE_VERSION} PARENT_SCOPE)
set(HYPRE_IS_PARALLEL  ${HYPRE_IS_PARALLEL} PARENT_SCOPE)
set(HYPRE_LIBRARY      ${HYPRE_LIBRARY} PARENT_SCOPE)
set(HYPRE_LIBRARIES    ${HYPRE_LIBRARIES} PARENT_SCOPE)
set(HYPRE_INCLUDE_DIR  ${HYPRE_INCLUDE_DIR} PARENT_SCOPE)
set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIRS} PARENT_SCOPE)

# HDF5
set(HDF5_FOUND          ${HDF5_FOUND} PARENT_SCOPE)
set(HDF5_BUILD_TARGET   ${HDF5_BUILD_TARGET} PARENT_SCOPE)
set(HDF5_VERSION        ${HDF5_VERSION} PARENT_SCOPE)
set(HDF5_C_LIBRARY      ${HDF5_C_LIBRARY} PARENT_SCOPE)
set(HDF5_HL_LIBRARY     ${HDF5_HL_LIBRARY} PARENT_SCOPE)
set(HDF5_C_LIBRARIES    ${HDF5_C_LIBRARIES} PARENT_SCOPE)
set(HDF5_LINK_LIBRARIES ${HDF5_LINK_LIBRARIES} PARENT_SCOPE)
set(HDF5_INCLUDE_DIR    ${HDF5_INCLUDE_DIR} PARENT_SCOPE)
set(HDF5_INCLUDE_DIRS   ${HDF5_INCLUDE_DIRS} PARENT_SCOPE)

# NETCDF
set(NETCDF_FOUND             ${NETCDF_FOUND} PARENT_SCOPE)
set(NETCDF_BUILD_TARGET      ${NETCDF_BUILD_TARGET} PARENT_SCOPE)
set(NETCDF_VERSION           ${NETCDF_VERSION} PARENT_SCOPE)
set(NETCDF_LARGE_DIMS        ${NETCDF_LARGE_DIMS} PARENT_SCOPE)
set(NETCDF_C_LIBRARY         ${NETCDF_C_LIBRARY} PARENT_SCOPE)
set(NETCDF_C_LIBRARIES       ${NETCDF_C_LIBRARIES} PARENT_SCOPE)
set(NETCDF_Fortran_LIBRARY   ${NETCDF_Fortran_LIBRARY} PARENT_SCOPE)
set(NETCDF_Fortran_LIBRARIES ${NETCDF_Fortran_LIBRARIES} PARENT_SCOPE)
set(NETCDF_INCLUDE_DIR       ${NETCDF_INCLUDE_DIR} PARENT_SCOPE)
set(NETCDF_INCLUDE_DIRS      ${NETCDF_INCLUDE_DIRS} PARENT_SCOPE)
