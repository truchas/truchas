# Set ENV PORTAGE_ROOT or PORTAGE_ROOT
# Defines imported library target portage

set(portage_search_path ${PORTAGE_ROOT} ENV PORTAGE_ROOT)
find_path(PORTAGE_INCLUDE_DIR NAMES portage/support/portage.h
    HINTS ${portage_search_path} PATH_SUFFIXES include)
find_library(PORTAGE_LIBRARY portage
    HINTS ${portage_search_path} PATH_SUFFIXES lib lib64)
find_library(WONTON_LIBRARY wonton
    HINTS ${portage_search_path} PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Portage DEFAULT_MSG
    PORTAGE_INCLUDE_DIR PORTAGE_LIBRARY WONTON_LIBRARY)

if(Portage_FOUND)
  set(PORTAGE_INCLUDE_DIRECTORIES "${PORTAGE_INCLUDE_DIR}")
  add_library(portage UNKNOWN IMPORTED)
  set_target_properties(portage PROPERTIES
      IMPORTED_LOCATION "${PORTAGE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${PORTAGE_INCLUDE_DIRECTORIES}")
  set_property(TARGET portage PROPERTY INTERFACE_LINK_LIBRARIES
     "${WONTON_LIBRARY}" "${MPI_CXX_LIBRARIES}")
endif()
