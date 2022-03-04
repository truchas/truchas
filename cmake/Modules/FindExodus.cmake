# This module finds the Exodus library and header file directory.
# EXODUS_FOUND is set to True if both are found and the following
# variables are returned.
#
#  EXODUS_INCLUDE_DIRS - required include directories
#  EXODUS_LIBRARIES    - required link libraries
#  EXODUS_VERSION      - version of exodus found
#
# This module also defines the imported library target "exodus".  It is
# generally enough to include "exodus" as a target link library; cmake
# will automatically handle adding the appropriate compile include flags
# and collection of link libraries.
#
# Set the variable CMAKE_PREFIX_PATH to provide a hint to the module for
# where to find the library and header file.  This is searched before the
# standard system locations.
#
# Find_package(NetCDF) should be run before using this module.

#if(Exodus_FIND_QUIETLY)
#  set(_find_netcdf_arg QUIET)
#endif()
#find_package(NetCDF ${_find_netcdf_arg})
#unset(_find_netcdf_arg)

find_path(EXODUS_INCLUDE_DIR exodusII.h)
find_library(EXODUS_LIBRARY exodus)

# Identify the ExodusII API version (not exactly the library version)
if(EXODUS_INCLUDE_DIR)
  set(exodusii_h ${EXODUS_INCLUDE_DIR}/exodusII.h)
  include(SearchHeaderFile)
  search_header_file(${exodusii_h} "EX_API_VERS" EXODUS_VERSION)
endif()

# Set EXODUS_FOUND to TRUE of the REQUIRED variables have been set.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Exodus
                                  REQUIRED_VARS EXODUS_LIBRARY EXODUS_INCLUDE_DIR
                                  VERSION_VAR EXODUS_VERSION)

if(EXODUS_FOUND)
  set(EXODUS_INCLUDE_DIRS ${EXODUS_INCLUDE_DIR} ${NETCDF_C_INCLUDE_DIRS})
  set(EXODUS_LIBRARIES ${EXODUS_LIBRARY} ${NETCDF_C_LIBRARIES})
  list(REMOVE_DUPLICATES EXODUS_INCLUDE_DIRS)
  mark_as_advanced(EXODUS_INCLUDE_DIR EXODUS_LIBRARY)
  if(NOT TARGET exodus)
    add_library(exodus UNKNOWN IMPORTED)
    set_target_properties(exodus PROPERTIES
        IMPORTED_LOCATION "${EXODUS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${EXODUS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${NETCDF_C_LIBRARIES}")
  endif()
endif()
