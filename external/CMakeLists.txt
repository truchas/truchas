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

# --- External projects built for Truchas

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

# GNU configure static/shared
set(TruchasExternal_SHARED_SWITCH --disable-shared)
if(BUILD_SHARED_LIBS)
  set(TruchasExternal_SHARED_SWITCH --enable-shared)
endif()

# Log all activity
set(TruchasExternal_LOG_OPTS
    LOG_DOWNLOAD  1
    LOG_UPDATE    1
    LOG_CONFIGURE 1
    LOG_BUILD     1
    LOG_TEST      1
    LOG_INSTALL   1
)



# HYPRE 
set(HYPRE_REQUIRED_VERSION 2.6.0b)
set(HYPRE_INSTALL_PREFIX ${TruchasExternal_INSTALL_PREFIX})
find_package(HYPRE)
if ( HYPRE_FOUND )
  if ( NOT ${HYPRE_VERSION} VERSION_EQUAL ${HYPRE_REQUIRED_VERSION} )
    message(STATUS "Incorrect HYPRE version ${HYPRE_VERSION} != ${HYPRE_REQUIRED_VERSION}")
    include(ExternalProject_HYPRE)
  elseif( ENABLE_MPI AND NOT HYPRE_IS_PARALLEL )
    message(STATUS "Require parallel HYPRE will add to the build")
    include(ExternalProject_HYPRE)
  else()  
    message(STATUS "Located HYPRE will not add to the build.")
    set(HYPRE_FOUND ${HYPRE_FOUND} PARENT_SCOPE)
    set(HYPRE_VERSION ${HYPRE_VERSION} PARENT_SCOPE)
    set(HYPRE_IS_PARALLEL ${HYPRE_IS_PARALLEL} PARENT_SCOPE)
    set(HYPRE_INCLUDE_DIR ${HYPRE_INCLUDE_DIR} PARENT_SCOPE)
    set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIRS} PARENT_SCOPE)
    set(HYPRE_LIBRARY ${HYPRE_LIBRARY} PARENT_SCOPE)
    set(HYPRE_LIBRARIES ${HYPRE_LIBRARIES} PARENT_SCOPE)
  endif()  
else()  
  message(STATUS "Could not locate HYPRE will add HYPRE to build.")
  include(ExternalProject_HYPRE)
endif()  

# HDF5
include(ExternalProject_ZLIB)
include(ExternalProject_HDF5)

# SWIG
find_package(SWIG)
if ( SWIG_FOUND )
  if  ( "${SWIG_VERSION}" VERSION_GREATER 2.0  )
    message(STATUS "Located Swig will not build.")
    #set(SWIG_FOUND ${SWIG_FOUND} PARENT_SCOPE)
    #set(SWIG_EXECUTABLE ${SWIG_EXECUTABLE} PARENT_SCOPE)
    #set(SWIG_VERSION ${SWIG_VERSION} PARENT_SCOPE)
  else()  
    message(STATUS "Swig version ${SWIG_VERSION} is not sufficient")
    include(ExternlProject_SWIG)
  endif()  
else()
  message(STATUS "Can not locate Swig will add Swig to build.")
  include(ExternalProject_SWIG)
endif()

# NetCDF
include(ExternalProject_NetCDF)




