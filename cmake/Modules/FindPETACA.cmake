# - Find PETACA
# ############################################################################ #
# Find the PETACA includes and library.
# Once done this will define
#
#  PETACA_MODULE_DIR     (PATH) - where to find the Fortran module files
#  PETACA_INCLUDE_DIRS   (PATH) - list of paths to include to compile
#  PETACA_LIBRARY        (FILE) - Petaca library
#  PETACA_LIBRARIES      (FILE) - list of libraries required to link against
#                                 Petaca
#  PETACA_FOUND          (BOOL) - True if Petaca found.
#
#  PETACA_VERSION        (STRING) - The version of Petaca found (x.y.z)
#
# This module defines the imported library target PETACA::PETACA and keeps the
# legacy imported target petaca for compatibility.
#
# Set PETACA_ROOT, PETACA_INSTALL_PREFIX, or CMAKE_PREFIX_PATH to provide a
# hint to the fallback search. The environment variable PETACA_ROOT is also
# honored.
#
# ############################################################################ #

find_package(PETACA CONFIG QUIET)
if(NOT PETACA_FOUND)
  find_package(Petaca CONFIG QUIET)
endif()

if(PETACA_FOUND OR Petaca_FOUND)
  set(_PETACA_CONFIG_TARGETS PETACA::PETACA Petaca::Petaca petaca::petaca petaca)
  foreach(_target IN LISTS _PETACA_CONFIG_TARGETS)
    if(TARGET ${_target})
      set(_PETACA_CONFIG_TARGET ${_target})
      break()
    endif()
  endforeach()

  if(_PETACA_CONFIG_TARGET)
    if(NOT TARGET PETACA::PETACA)
      add_library(PETACA::PETACA INTERFACE IMPORTED)
      target_link_libraries(PETACA::PETACA INTERFACE ${_PETACA_CONFIG_TARGET})
    endif()

    if(NOT TARGET petaca)
      add_library(petaca INTERFACE IMPORTED)
      target_link_libraries(petaca INTERFACE PETACA::PETACA)
    endif()

    get_target_property(_PETACA_IMPORTED_LOCATION PETACA::PETACA IMPORTED_LOCATION)
    get_target_property(_PETACA_INCLUDE_DIRS PETACA::PETACA INTERFACE_INCLUDE_DIRECTORIES)
    if(_PETACA_IMPORTED_LOCATION)
      set(PETACA_LIBRARY ${_PETACA_IMPORTED_LOCATION})
    endif()
    if(_PETACA_INCLUDE_DIRS)
      set(PETACA_MODULE_DIR ${_PETACA_INCLUDE_DIRS})
      set(PETACA_INCLUDE_DIRS ${_PETACA_INCLUDE_DIRS})
    endif()
    set(PETACA_LIBRARIES PETACA::PETACA)
    if(NOT PETACA_VERSION AND Petaca_VERSION)
      set(PETACA_VERSION ${Petaca_VERSION})
    endif()
    set(PETACA_FOUND True)

    unset(_PETACA_CONFIG_TARGET)
    unset(_PETACA_CONFIG_TARGETS)
    unset(_PETACA_IMPORTED_LOCATION)
    unset(_PETACA_INCLUDE_DIRS)
    return()
  endif()

  unset(_PETACA_CONFIG_TARGETS)
endif()

if(PETACA_FIND_QUIETLY)
  set(_FIND_YAJL_ARG QUIET)
endif()
find_package(YAJL "2.0.4" ${_FIND_YAJL_ARG})

set(_PETACA_ROOTS ${PETACA_ROOT} ${PETACA_INSTALL_PREFIX} ENV PETACA_ROOT)

if(YAJL_FOUND)
  find_library(PETACA_LIBRARY petaca
    HINTS ${_PETACA_ROOTS}
    PATH_SUFFIXES lib lib64)

  # Module files are installed with the library file.
  if(PETACA_LIBRARY)
    get_filename_component(PETACA_MODULE_DIR ${PETACA_LIBRARY} DIRECTORY)
  endif()

  # No version number is currently available
  set(PETACA_VERSION PETACA_VERSION-NOTFOUND)

  if(PETACA_LIBRARY AND PETACA_MODULE_DIR)
    set(PETACA_LIBRARIES ${PETACA_LIBRARY} ${YAJL_LIBRARIES})
    set(PETACA_INCLUDE_DIRS ${PETACA_MODULE_DIR})

    if(NOT TARGET PETACA::PETACA)
      add_library(PETACA::PETACA UNKNOWN IMPORTED)
      set_target_properties(PETACA::PETACA PROPERTIES
          IMPORTED_LOCATION "${PETACA_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${PETACA_INCLUDE_DIRS}")
      target_link_libraries(PETACA::PETACA INTERFACE YAJL::YAJL)
    endif()

    if(NOT TARGET petaca)
      add_library(petaca INTERFACE IMPORTED)
      target_link_libraries(petaca INTERFACE PETACA::PETACA)
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETACA
  REQUIRED_VARS PETACA_LIBRARY PETACA_MODULE_DIR YAJL_FOUND
  VERSION_VAR PETACA_VERSION
  REASON_FAILURE_MESSAGE "PETACA requires YAJL version 2.0.4 or newer.")

mark_as_advanced(PETACA_LIBRARY PETACA_MODULE_DIR)
unset(_FIND_YAJL_ARG)
unset(_PETACA_ROOTS)
