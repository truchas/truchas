# Find MUMPS.
#
# Prefer an installed MUMPS CMake package when one is available, but support
# direct MUMPS installs that only provide headers and libraries.

set(_MUMPS_config_args)
if(MUMPS_FIND_VERSION)
  list(APPEND _MUMPS_config_args "${MUMPS_FIND_VERSION}")
endif()
list(APPEND _MUMPS_config_args CONFIG QUIET)
if(MUMPS_FIND_COMPONENTS)
  list(APPEND _MUMPS_config_args COMPONENTS ${MUMPS_FIND_COMPONENTS})
endif()
find_package(MUMPS ${_MUMPS_config_args})
unset(_MUMPS_config_args)

if(MUMPS_FOUND)
  if(NOT MUMPS_VERSION AND MUMPS_UPSTREAM_VERSION)
    set(MUMPS_VERSION "${MUMPS_UPSTREAM_VERSION}")
  endif()
  return()
endif()

set(_MUMPS_search_paths
  ${MUMPS_ROOT}
  ${MUMPS_INSTALL_PREFIX}
  ${MUMPS_DIR}
  ENV MUMPS_ROOT
)

if(MUMPS_FIND_COMPONENTS)
  set(_MUMPS_components ${MUMPS_FIND_COMPONENTS})
else()
  set(_MUMPS_components d z)
endif()

find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_struc.h
  HINTS ${_MUMPS_search_paths}
  PATH_SUFFIXES include)

find_library(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  HINTS ${_MUMPS_search_paths}
  PATH_SUFFIXES lib lib64)
find_library(MUMPS_PORD_LIBRARY
  NAMES pord
  HINTS ${_MUMPS_search_paths}
  PATH_SUFFIXES lib lib64)
find_library(MUMPS_SCALAPACK_LIBRARY
  NAMES scalapack
  HINTS ${_MUMPS_search_paths}
  PATH_SUFFIXES lib lib64)

foreach(_MUMPS_comp IN LISTS _MUMPS_components)
  string(TOLOWER "${_MUMPS_comp}" _MUMPS_comp_lc)
  string(TOUPPER "${_MUMPS_comp}" _MUMPS_comp_uc)
  find_library(MUMPS_${_MUMPS_comp_uc}_LIBRARY
    NAMES ${_MUMPS_comp_lc}mumps
    HINTS ${_MUMPS_search_paths}
    PATH_SUFFIXES lib lib64)
  if(MUMPS_${_MUMPS_comp_uc}_LIBRARY)
    set(MUMPS_${_MUMPS_comp_lc}_FOUND TRUE)
  else()
    set(MUMPS_${_MUMPS_comp_lc}_FOUND FALSE)
  endif()
endforeach()

if(MUMPS_INCLUDE_DIR)
  include(SearchHeaderFile)
  search_header_file(${MUMPS_INCLUDE_DIR}/dmumps_struc.h
    MUMPS_RELEASE_VERSION MUMPS_VERSION STRIP_QUOTES)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_COMMON_LIBRARY MUMPS_INCLUDE_DIR
  VERSION_VAR MUMPS_VERSION
  HANDLE_COMPONENTS)

if(MUMPS_FOUND)
  set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
  set(MUMPS_LIBRARIES ${MUMPS_COMMON_LIBRARY})

  if(MUMPS_PORD_LIBRARY AND NOT TARGET MUMPS::PORD)
    add_library(MUMPS::PORD UNKNOWN IMPORTED)
    set_target_properties(MUMPS::PORD PROPERTIES
      IMPORTED_LOCATION "${MUMPS_PORD_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
  endif()

  if(NOT TARGET MUMPS::COMMON)
    add_library(MUMPS::COMMON UNKNOWN IMPORTED)
    set_target_properties(MUMPS::COMMON PROPERTIES
      IMPORTED_LOCATION "${MUMPS_COMMON_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
    if(MUMPS_PORD_LIBRARY)
      target_link_libraries(MUMPS::COMMON INTERFACE MUMPS::PORD)
      list(APPEND MUMPS_LIBRARIES ${MUMPS_PORD_LIBRARY})
    endif()
    if(MUMPS_SCALAPACK_LIBRARY)
      target_link_libraries(MUMPS::COMMON INTERFACE "${MUMPS_SCALAPACK_LIBRARY}")
      list(APPEND MUMPS_LIBRARIES ${MUMPS_SCALAPACK_LIBRARY})
    endif()
    target_link_libraries(MUMPS::COMMON INTERFACE LAPACK::LAPACK MPI::MPI_Fortran)
  endif()

  foreach(_MUMPS_comp IN LISTS _MUMPS_components)
    string(TOLOWER "${_MUMPS_comp}" _MUMPS_comp_lc)
    string(TOUPPER "${_MUMPS_comp}" _MUMPS_comp_uc)
    if(MUMPS_${_MUMPS_comp_uc}_LIBRARY AND NOT TARGET MUMPS::${_MUMPS_comp_uc}MUMPS)
      add_library(MUMPS::${_MUMPS_comp_uc}MUMPS UNKNOWN IMPORTED)
      set_target_properties(MUMPS::${_MUMPS_comp_uc}MUMPS PROPERTIES
        IMPORTED_LOCATION "${MUMPS_${_MUMPS_comp_uc}_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES MUMPS::COMMON)
      list(APPEND MUMPS_LIBRARIES ${MUMPS_${_MUMPS_comp_uc}_LIBRARY})
    endif()
  endforeach()

  if(NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS INTERFACE IMPORTED)
    foreach(_MUMPS_comp IN LISTS _MUMPS_components)
      string(TOUPPER "${_MUMPS_comp}" _MUMPS_comp_uc)
      if(TARGET MUMPS::${_MUMPS_comp_uc}MUMPS)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::${_MUMPS_comp_uc}MUMPS)
      endif()
    endforeach()
  endif()
endif()

mark_as_advanced(
  MUMPS_INCLUDE_DIR
  MUMPS_COMMON_LIBRARY
  MUMPS_PORD_LIBRARY
  MUMPS_SCALAPACK_LIBRARY
)
foreach(_MUMPS_comp IN LISTS _MUMPS_components)
  string(TOUPPER "${_MUMPS_comp}" _MUMPS_comp_uc)
  mark_as_advanced(MUMPS_${_MUMPS_comp_uc}_LIBRARY)
endforeach()
unset(_MUMPS_comp)
unset(_MUMPS_comp_lc)
unset(_MUMPS_comp_uc)
unset(_MUMPS_components)
unset(_MUMPS_search_paths)
