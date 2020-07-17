# Defines imported library target chaco
# Set environment variable CHACO_ROOT or CMake variable CHACO_ROOT

set(chaco_search_path ${CHACO_ROOT} ENV CHACO_ROOT)
find_library(CHACO_LIBRARY chaco
    HINTS ${chaco_search_path} PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Chaco DEFAULT_MSG CHACO_LIBRARY)

if(Chaco_FOUND)
  add_library(chaco UNKNOWN IMPORTED)
  set_target_properties(chaco PROPERTIES
      IMPORTED_LOCATION "${CHACO_LIBRARY}")
  set_property(TARGET chaco PROPERTY INTERFACE_LINK_LIBRARIES
     "${CHACO_LIBRARY}")
endif()
