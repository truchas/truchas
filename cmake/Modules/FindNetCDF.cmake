# - Find NetCDF
# Find the NetCDF includes and libraries
#
# This module finds the netCDF library and header file directory for the
# C interface, and if the "Fortran" component is requested, those for the
# Fortran interface.
#
# If everything is found, NETCDF_FOUND is set to True and the following
# variables are defined.
#
#  NETCDF_VERSION              - Version of NetCDF found
#  NETCDF_HAS_NC4              - True if NetCDF includes NetCDF-4 support
#  NETCDF_LARGE_MODEL          - True if header file has large model tweaks
#  NETCDF_LIBRARIES            - Link libraries for all interfaces
#  NETCDF_INCLUDE_DIRS         - All include directories
#  NETCDF_C_LIBRARIES          - Link libraries for C interface
#  NETCDF_C_INCLUDE_DIRS       - Include directories for C interface
#
# If the Fortran interface was requested, these are also defined.
#
#  NETCDF_Fortran_LIBRARIES    - Link libraries for Fortran interface
#  NETCDF_Fortran_INCLUDE_DIRS - Include directories for Fortran interface
#
# The NETCDF_LIBARIES and NETCDF_INCLUDE_DIRS variables are just unions of the
# correspondig language-specific variables.

# Validate the optional component list
set(NETCDF_VALID_COMPONENTS Fortran)
set(NETCDF_LANGUAGE_BINDINGS "C")
if(NetCDF_FIND_COMPONENTS)
  foreach(_comp ${NetCDF_FIND_COMPONENTS})
    list(FIND NETCDF_VALID_COMPONENTS ${_comp} _loc)
    if(${_loc} EQUAL -1)
      message(FATAL_ERROR "\"${_comp}\" is not a valid NetCDF component.")
    else()
      list(APPEND NETCDF_LANGUAGE_BINDINGS ${_comp})
    endif()
  endforeach()
  unset(_comp)
  unset(_loc)
endif()
unset(NETCDF_VALID_COMPONENTS)

# Check whether netCDF-4 support is available (uses HDF5)
find_program(NC_CONFIG nc-config)
if(NC_CONFIG)
  execute_process(COMMAND ${NC_CONFIG} --has-nc4
                  OUTPUT_VARIABLE NETCDF_HAS_NC4
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(NETCDF_HAS_NC4 NETCDF_HAS_NC4-NOTFOUND)
endif()

# Should call find_package(HDF5) before calling this module
#if(NETCDF_HAS_NC4)
#  find_package(HDF5)
#endif()

# Always find the C interface
find_path(NETCDF_C_INCLUDE_DIR netcdf.h HINTS ENV NETCDF_DIR PATH_SUFFIXES include)
find_library(NETCDF_C_LIBRARY netcdf HINTS ENV NETCDF_DIR PATH_SUFFIXES lib)
mark_as_advanced(NETCDF_C_INCLUDE_DIR NETCDF_C_LIBRARY)
set(NETCDF_C_INCLUDE_DIRS ${NETCDF_C_INCLUDE_DIR})
set(NETCDF_C_LIBRARIES ${NETCDF_C_LIBRARY})
set(NETCDF_INCLUDE_DIRS ${NETCDF_C_INCLUDE_DIRS})
set(NETCDF_LIBRARIES ${NETCDF_C_LIBRARIES})

# Fortran interface (optional)
list(FIND NETCDF_LANGUAGE_BINDINGS Fortran index)
if(index GREATER -1)
  find_path(NETCDF_Fortran_INCLUDE_DIR NAMES netcdf.mod NETCDF.mod
    HINTS ENV NETCDF_DIR PATH_SUFFIXES include)
  find_library(NETCDF_Fortran_LIBRARY netcdff HINTS ENV NETCDF_DIR PATH_SUFFIXES lib)
  mark_as_advanced(NETCDF_Fortran_INCLUDE_DIR NETCDF_Fortran_LIBRARY)
  if(NETCDF_Fortran_INCLUDE_DIR AND NETCDF_Fortran_LIBRARY)
    set(NetCDF_Fortran_FOUND True)
  else()
    set(NetCDF_Fortran_FOUND False)
  endif()
  set(NETCDF_Fortran_INCLUDE_DIRS ${NETCDF_Fortran_INCLUDE_DIR})
  set(NETCDF_Fortran_LIBRARIES ${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARIES})
  list(APPEND NETCDF_INCLUDE_DIRS ${NETCDF_Fortran_INCLUDE_DIR})
  list(REMOVE_DUPLICATES NETCDF_INCLUDE_DIRS)
  list(INSERT NETCDF_LIBRARIES 0 ${NETCDF_Fortran_LIBRARY})
endif()

# NETCDF Version
if(NC_CONFIG)
  execute_process(COMMAND ${NC_CONFIG} --version
                  OUTPUT_VARIABLE NETCDF_VERSION)
  string(REGEX REPLACE "[\n\r ]" "" NETCDF_VERSION "${NETCDF_VERSION}")
  string(REGEX REPLACE "netCDF"  "" NETCDF_VERSION "${NETCDF_VERSION}")
else()
  set(NETCDF_VERSION NETCDF_VERSION-NOTFOUND)
endif()

# Check for the "large model modifications" desired by the ExodusII library
if(NETCDF_C_INCLUDE_DIR)
  set(netcdf_h "${NETCDF_C_INCLUDE_DIR}/netcdf.h")
  include(SearchHeaderFile)
  search_header_file(${netcdf_h} "NC_MAX_DIMS" max_dims)
  search_header_file(${netcdf_h} "NC_MAX_VARS" max_vars)
  search_header_file(${netcdf_h} "NC_MAX_VAR_DIMS" max_var_dims)
  if((max_dims EQUAL 65536) AND (max_vars EQUAL 524288) AND (max_var_dims EQUAL 8))
    set(NETCDF_LARGE_MODEL True)
  else()
    set(NETCDF_LARGE_MODEL False)
  endif()
  unset(max_dims)
  unset(max_vars)
  unset(max_var_dims)
  unset(netcdf_h)
else()
  set(NETCDF_LARGE_MODEL NETCDF_LARGE_MODEL-NOTFOUND)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS
    VERSION_VAR NETCDF_VERSION
    HANDLE_COMPONENTS)

if(NETCDF_FOUND)
  if(NETCDF_C_LIBRARY AND NETCDF_C_INCLUDE_DIR)
    if(NOT TARGET netcdf)
      add_library(netcdf UNKNOWN IMPORTED)
      set_target_properties(netcdf PROPERTIES
          IMPORTED_LOCATION "${NETCDF_C_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_C_INCLUDE_DIR}")
      if(NETCDF_HAS_NC4)
        set_target_properties(netcdf PROPERTIES
            INTERFACE_LINK_LIBRARIES "${HDF5_C_LIBRARIES}")
      endif()
    endif()
  endif()
  if(NETCDF_Fortran_LIBRARY AND NETCDF_Fortran_INCLUDE_DIR)
    if(NOT TARGET netcdff)
      add_library(netcdff UNKNOWN IMPORTED)
      set_target_properties(netcdff PROPERTIES
          IMPORTED_LOCATION "${NETCDF_Fortran_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_Fortran_INCLUDE_DIR}"
          INTERFACE_LINK_LIBRARIES netcdf)
    endif()
  endif()
endif()

# clean up
unset(NETCDF_LANGUAGE_BINDINGS)
unset(NC_CONFIG)
