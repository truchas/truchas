# - Find the HYPRE (High Performance Preconditioners) Implementation
# HYPRE is a popular software package that implements several scalable
# linear solvers. This module searches for a HYPRE installation and
# returns the pointers to the include path and library
#
# Defines imported library target chaparral
# Set environment variable CHAPARRAL_ROOT or CMake variable CHAPARRAL_ROOT

set(chaparral_search_paths ${CHAPARRAL_ROOT} ENV CHAPARRAL_ROOT)
find_library(CHAPARRAL_LIBRARY VF
    HINTS ${chaparral_search_paths} PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Chaparral DEFAULT_MSG CHAPARRAL_LIBRARY)

if(Chaparral_FOUND)
  set(CHAPARRAL_LIBRARIES ${CHAPARRAL_LIBRARY} ${MPI_C_LIBRARIES})
  if(NOT TARGET chaparral)
    add_library(chaparral UNKNOWN IMPORTED)
    set_target_properties(chaparral PROPERTIES
        IMPORTED_LOCATION "${CHAPARRAL_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "${CHAPARRAL_LIBRARIES}")
  endif()
endif()
