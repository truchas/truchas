# Defines imported library target scorpio
# Set environment variable SCORPIO_ROOT or CMake variable SCORPIO_ROOT

set(scorpio_search_paths ${SCORPIO_ROOT} ENV SCORPIO_ROOT)
find_library(SCORPIO_LIBRARY scorpio
    HINTS ${scorpio_search_paths} PATH_SUFFIXES lib lib64)
find_path(SCORPIO_INCLUDE_DIR NAMES scorpio.h
    HINTS ${scorpio_search_paths} PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Scorpio DEFAULT_MSG
    SCORPIO_LIBRARY SCORPIO_INCLUDE_DIR)

if(Scorpio_FOUND)
  set(SCORPIO_INCLUDE_DIRS ${SCORPIO_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})
  set(SCORPIO_LIBRARIES ${SCORPIO_LIBRARY})
  if(NOT TARGET scorpio)
    add_library(scorpio UNKNOWN IMPORTED)
    set_target_properties(scorpio PROPERTIES
        IMPORTED_LOCATION "${SCORPIO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SCORPIO_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${SCORPIO_LIBRARIES}")
  endif()
endif()
