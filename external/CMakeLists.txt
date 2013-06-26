# #  -*- mode: cmake -*-
# ############################################################################ #
#                                                                              #
# TruchasExternal                                                              #
#                                                                              #
# ############################################################################ #

project(TruchasExternal)

# --- Add the local cmake directory 
list(APPEND CMAKE_MODULE_PATH ${TruchasExternal_SOURCE_DIR}/cmake)

# --- Define the external project installation  
if (NOT TruchasExternal_INSTALL_PREFIX)
  set(TruchasExternal_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/external
      CACHE PATH "Truchas external projects install prefix")
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
#BROKENfind_package(HYPRE)

#BROKENif ( NOT HYPRE_FOUND )
#BROKEN  include(ExternalProject_HYPRE)
#BORKENendif()  
include(ExternalProject_HYPRE)

# HDF5
include(ExternalProject_ZLIB)
include(ExternalProject_HDF5)

# SWIG
include(ExternalProject_SWIG)

# NetCDF
include(ExternalProject_NetCDF)




