# Defines imported library target Scorpio::Scorpio
# Set environment variable SCORPIO_ROOT or CMake variable SCORPIO_ROOT

find_package(Scorpio CONFIG QUIET)

if(Scorpio_FOUND)
  if(TARGET Scorpio::Scorpio AND NOT TARGET scorpio)
    add_library(scorpio INTERFACE IMPORTED)
    target_link_libraries(scorpio INTERFACE Scorpio::Scorpio)
  endif()
  return()
endif()

set(_Scorpio_ROOTS ${Scorpio_ROOT} ${SCORPIO_ROOT} ENV Scorpio_ROOT ENV SCORPIO_ROOT)
find_library(SCORPIO_LIBRARY scorpio
    HINTS ${_Scorpio_ROOTS} PATH_SUFFIXES lib lib64)
find_path(SCORPIO_INCLUDE_DIR NAMES scorpio.h
    HINTS ${_Scorpio_ROOTS} PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Scorpio
    REQUIRED_VARS SCORPIO_LIBRARY SCORPIO_INCLUDE_DIR)

if(Scorpio_FOUND)
  set(SCORPIO_INCLUDE_DIRS ${SCORPIO_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
  set(SCORPIO_LIBRARIES ${SCORPIO_LIBRARY})

  if(NOT TARGET Scorpio::Scorpio)
    add_library(Scorpio::Scorpio UNKNOWN IMPORTED)
    set_target_properties(Scorpio::Scorpio PROPERTIES
        IMPORTED_LOCATION "${SCORPIO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SCORPIO_INCLUDE_DIR}")
    if(TARGET MPI::MPI_C)
      target_link_libraries(Scorpio::Scorpio INTERFACE MPI::MPI_C)
    endif()
  endif()

  if(NOT TARGET scorpio)
    add_library(scorpio INTERFACE IMPORTED)
    target_link_libraries(scorpio INTERFACE Scorpio::Scorpio)
  endif()
endif()

mark_as_advanced(SCORPIO_INCLUDE_DIR SCORPIO_LIBRARY)
unset(_Scorpio_ROOTS)
