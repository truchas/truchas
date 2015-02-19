#
# Find the Exodus library and include directory.
# Variables set by this find module:
#
#  EXODUS_FOUND        - true if exodus was found
#  EXODUS_INCLUDE_DIRS - where to find exodusII.h
#  EXODUS_LIBRARIES    - list of link libraries for exodus
#  EXODUS_VERSION      - version of exodus found
#
# This module will first search the Exodus installation root specified
# by the environment variable EXODUS_ROOT or the cmake variable 
# EXODUS_INSTALL_PREFIX for the library and header file before searching
# in the standards locations.
#

# Search paths
set(_exodus_search_paths
    ${EXODUS_INSTALL_PREFIX}
    ENV EXODUS_ROOT)

# Locate the include directory
find_path(EXODUS_INCLUDE_DIR
          NAMES exodusII.h
          HINTS ${_exodus_search_paths}
          PATH_SUFFIXES include)

# Locate the library        
find_library(EXODUS_LIBRARY
             NAMES exoIIv2c
             HINTS ${_exodus_search_paths}
             PATH_SUFFIXES lib)

# Identify the library version
find_file(_exodus_h NAMES exodusII.h HINTS ${EXODUS_INCLUDE_DIR})
if(_exodus_h)
  file(STRINGS "${_exodus_h}" _ex_api_vers_nodot REGEX "^#define EX_API_VERS_NODOT")
  string(REGEX REPLACE "[^0-9]" "" EXODUS_VERSION "${_ex_api_vers_nodot}")
else()
  unset(EXODUS_VERSION)
endif()

# Set EXODUS_FOUND to TRUE of the REQUIRED variables have been set.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Exodus
                                  REQUIRED_VARS EXODUS_LIBRARY EXODUS_INCLUDE_DIR 
                                  VERSION_VAR EXODUS_VERSION
                                  HANDLE_COMPONENTS)
if(EXODUS_FOUND)
  set(EXODUS_INCLUDE_DIRS ${EXODUS_INCLUDE_DIR})
  set(EXODUS_LIBRARIES ${EXODUS_LIBRARY})
  if(NETCDF_FOUND AND NETCDF_C_LIBRARIES)
    list(APPEND EXODUS_LIBRARIES ${NETCDF_C_LIBRARIES})
  else()
    message(WARNING "EXODUS_LIBRARIES does not contain the NetCDF libraries")
  endif()
endif()
