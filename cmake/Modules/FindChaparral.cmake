# Defines imported library target chaparral
# Set environment variable CHAPARRAL_ROOT or CMake variable CHAPARRAL_ROOT

set(chaparral_search_paths ${CHAPARRAL_ROOT} ENV CHAPARRAL_ROOT)
find_library(CHAPARRAL_LIBRARY VF
    HINTS ${chaparral_search_paths} PATH_SUFFIXES lib lib64)
find_path(CHAPARRAL_INCLUDE_DIR NAMES vf_api.h
    HINTS ${chaparral_search_paths} PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Chaparral DEFAULT_MSG
    CHAPARRAL_LIBRARY CHAPARRAL_INCLUDE_DIR)

if(Chaparral_FOUND)
  set(CHAPARRAL_INCLUDE_DIRS ${CHAPARRAL_INCLUDE_DIR})
  set(CHAPARRAL_LIBRARIES ${CHAPARRAL_LIBRARY} ${MPI_C_LIBRARIES})
  if(NOT TARGET chaparral)
    add_library(chaparral UNKNOWN IMPORTED)
    set_target_properties(chaparral PROPERTIES
        IMPORTED_LOCATION "${CHAPARRAL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CHAPARRAL_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${CHAPARRAL_LIBRARIES}")
  endif()
endif()
