# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindMETIS
-------
Michael Hirsch, Ph.D.

Finds the METIS library.
NOTE: If libparmetis used, libmetis must also be linked.

Imported Targets
^^^^^^^^^^^^^^^^

METIS::METIS

Result Variables
^^^^^^^^^^^^^^^^

METIS_LIBRARIES
  libraries to be linked

METIS_INCLUDE_DIRS
  dirs to be included

#]=======================================================================]


set(_METIS_search_paths
  ${METIS_ROOT}
  ${METIS_DIR}
  ENV METIS_ROOT
)

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  HINTS ${_METIS_search_paths}
  PATH_SUFFIXES include METIS openmpi-x86_64 mpich-x86_64
  DOC "METIS include directory")

find_library(METIS_LIBRARY
  NAMES metis
  HINTS ${_METIS_search_paths}
  PATH_SUFFIXES lib lib64 METIS libmetis
  DOC "METIS library")

if(METIS_INCLUDE_DIR)
  file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _METIS_version_lines
    REGEX "^#define METIS_VER_(MAJOR|MINOR|SUBMINOR)[ \t]+[0-9]+")
  foreach(_METIS_version_line IN LISTS _METIS_version_lines)
    if(_METIS_version_line MATCHES "^#define METIS_VER_MAJOR[ \t]+([0-9]+)")
      set(_METIS_version_major "${CMAKE_MATCH_1}")
    elseif(_METIS_version_line MATCHES "^#define METIS_VER_MINOR[ \t]+([0-9]+)")
      set(_METIS_version_minor "${CMAKE_MATCH_1}")
    elseif(_METIS_version_line MATCHES "^#define METIS_VER_SUBMINOR[ \t]+([0-9]+)")
      set(_METIS_version_subminor "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(DEFINED _METIS_version_major AND DEFINED _METIS_version_minor AND DEFINED _METIS_version_subminor)
    set(METIS_VERSION "${_METIS_version_major}.${_METIS_version_minor}.${_METIS_version_subminor}")
  endif()
endif()

if(ParMETIS IN_LIST METIS_FIND_COMPONENTS)
  find_path(PARMETIS_INCLUDE_DIR
    NAMES parmetis.h
    HINTS ${_METIS_search_paths}
    PATH_SUFFIXES include METIS openmpi-x86_64 mpich-x86_64
    DOC "ParMETIS include directory")
  find_library(PARMETIS_LIBRARY
    NAMES parmetis
    HINTS ${_METIS_search_paths}
    PATH_SUFFIXES lib lib64 METIS libmetis
    DOC "ParMETIS library")
  if(PARMETIS_INCLUDE_DIR AND PARMETIS_LIBRARY)
    set(METIS_ParMETIS_FOUND TRUE)
  else()
    set(METIS_ParMETIS_FOUND FALSE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  REQUIRED_VARS METIS_LIBRARY METIS_INCLUDE_DIR
  VERSION_VAR METIS_VERSION
  HANDLE_COMPONENTS)

if(METIS_FOUND)
  set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
  set(METIS_LIBRARIES ${METIS_LIBRARY})

  if(NOT TARGET METIS::METIS)
    add_library(METIS::METIS UNKNOWN IMPORTED)
    set_target_properties(METIS::METIS PROPERTIES
      IMPORTED_LOCATION "${METIS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")
  endif()

  if(METIS_ParMETIS_FOUND)
    list(APPEND METIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
    list(REMOVE_DUPLICATES METIS_INCLUDE_DIRS)
    list(PREPEND METIS_LIBRARIES ${PARMETIS_LIBRARY})
    if(NOT TARGET METIS::ParMETIS)
      add_library(METIS::ParMETIS UNKNOWN IMPORTED)
      set_target_properties(METIS::ParMETIS PROPERTIES
        IMPORTED_LOCATION "${PARMETIS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES METIS::METIS)
    endif()
  endif()

  message(VERBOSE "METIS libraries: ${METIS_LIBRARIES}
  METIS include directories: ${METIS_INCLUDE_DIRS}")
endif()

mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY)
unset(_METIS_search_paths)
unset(_METIS_version_line)
unset(_METIS_version_lines)
unset(_METIS_version_major)
unset(_METIS_version_minor)
unset(_METIS_version_subminor)
