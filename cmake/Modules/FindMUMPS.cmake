# Call this module with FIND_PACKAGE(MUMPS) in any CMakeLists.txt
# that requires MUMPS include directories or libraries. The search
# can be bypassed by explicitly setting MUMPS_INCLUDE_DIRS and
# MUMPS_LIBRARIES variables.
#
# The variables MUMPS_INSTALL_PREFIX or ENV MUMPS_ROOT control
# search locations for include and libraries.

find_package(MUMPS CONFIG)

if(MUMPS_FOUND)
  message("Found MUMPS CMake config")
else()
  message("Did not find MUMPS CMake config. Manually configuring.")
  set(mumps_search_paths ${MUMPS_DIR} ENV MUMPS_ROOT)
  find_library(MUMPS_LIBRARY mumps_common HINTS ${mumps_search_paths} PATH_SUFFIXES lib)
  find_library(DMUMPS_LIBRARY dmumps HINTS ${mumps_search_paths} PATH_SUFFIXES lib)
  find_library(ZMUMPS_LIBRARY zmumps HINTS ${mumps_search_paths} PATH_SUFFIXES lib)
  find_library(SCALAPACK_LIBRARY scalapack HINTS ${mumps_search_paths} PATH_SUFFIXES lib)
  find_library(PORD_LIBRARY pord HINTS ${mumps_search_paths} PATH_SUFFIXES lib)
  find_path(MUMPS_INCLUDE_DIR dmumps_struc.h HINTS ${mumps_search_paths} PATH_SUFFIXES include)

  if(MUMPS_INCLUDE_DIR)
    include(SearchHeaderFile)
    search_header_file(${MUMPS_INCLUDE_DIR}/dmumps_struc.h
      MUMPS_RELEASE_VERSION MUMPS_VERSION STRIP_QUOTES)
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MUMPS
    REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR
    VERSION_VAR MUMPS_VERSION)

  if(MUMPS_FOUND)
    set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
    set(MUMPS_LIBRARIES ${MUMPS_LIBRARY} ${DMUMPS_LIBRARY} ${ZMUMPS_LIBRARY} ${SCALAPACK_LIBRARY} ${PORD_LIBRARY})

    if(NOT TARGET MUMPS::DMUMPS)
      add_library(MUMPS::DMUMPS UNKNOWN IMPORTED)
      set_target_properties(MUMPS::DMUMPS PROPERTIES
        IMPORTED_LOCATION "${DMUMPS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${DMUMPS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${DMUMPS_LIBRARIES};${MUMPS_LIBRARIES};${SCALAPACK_LIBRARIES};${PORD_LIBRARIES}")
    endif()

    if(NOT TARGET MUMPS::ZMUMPS)
      add_library(MUMPS::ZMUMPS UNKNOWN IMPORTED)
      set_target_properties(MUMPS::ZMUMPS PROPERTIES
        IMPORTED_LOCATION "${ZMUMPS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ZMUMPS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${ZMUMPS_LIBRARIES}")
    endif()
  endif()
endif()
