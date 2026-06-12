# Defines imported library target Chaparral::Chaparral
# Set environment variable CHAPARRAL_ROOT or CMake variable CHAPARRAL_ROOT

find_package(Chaparral CONFIG QUIET)

if(Chaparral_FOUND)
  if(TARGET Chaparral::Chaparral AND NOT TARGET chaparral)
    add_library(chaparral INTERFACE IMPORTED)
    target_link_libraries(chaparral INTERFACE Chaparral::Chaparral)
  endif()
  return()
endif()

set(_Chaparral_ROOTS ${Chaparral_ROOT} ${CHAPARRAL_ROOT} ENV Chaparral_ROOT ENV CHAPARRAL_ROOT)
find_library(CHAPARRAL_LIBRARY VF
    HINTS ${_Chaparral_ROOTS} PATH_SUFFIXES lib lib64)
find_path(CHAPARRAL_INCLUDE_DIR NAMES vf_api.h
    HINTS ${_Chaparral_ROOTS} PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Chaparral
    REQUIRED_VARS CHAPARRAL_LIBRARY CHAPARRAL_INCLUDE_DIR)

if(Chaparral_FOUND)
  set(CHAPARRAL_INCLUDE_DIRS ${CHAPARRAL_INCLUDE_DIR})
  set(CHAPARRAL_LIBRARIES ${CHAPARRAL_LIBRARY} ${MPI_C_LIBRARIES})

  if(NOT TARGET Chaparral::Chaparral)
    add_library(Chaparral::Chaparral UNKNOWN IMPORTED)
    set_target_properties(Chaparral::Chaparral PROPERTIES
        IMPORTED_LOCATION "${CHAPARRAL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CHAPARRAL_INCLUDE_DIRS}")
    if(TARGET MPI::MPI_C)
      target_link_libraries(Chaparral::Chaparral INTERFACE MPI::MPI_C)
    endif()
  endif()

  if(NOT TARGET chaparral)
    add_library(chaparral INTERFACE IMPORTED)
    target_link_libraries(chaparral INTERFACE Chaparral::Chaparral)
  endif()
endif()

mark_as_advanced(CHAPARRAL_INCLUDE_DIR CHAPARRAL_LIBRARY)
unset(_Chaparral_ROOTS)
