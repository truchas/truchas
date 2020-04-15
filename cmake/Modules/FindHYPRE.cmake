# - Find the HYPRE (High Performance Preconditioners) Implementation
# HYPRE is a popular software package that implements several scalable
# linear solvers. This module searches for a HYPRE installation and
# returns the pointers to the include path and library
#
# === Variables ===
#
# HYPRE_FOUND             Flag indicating if HYPRE is found
# HYPRE_VERSION           HYPRE version
# HYPRE_INCLUDE_DIRS      List of paths required to compile
#                         with HYPRE
# HYPRE_LIBRARIES         List libraries (targets?) required
#                         to link against HYPRE.
# HYPRE_IS_PARALLEL       HYPRE was found
#
# HYPRE_INSTALL_PREFIX    Hypre installation prefix. Used in the
#                         library and include file search.
# === Usage ===
#
# Call this module with FIND_PACKAGE(HYPRE) in any CMakeLists.txt
# that requires HYPRE include directories or libraries. The search
# can be bypassed by explicitly setting HYPRE_INCLUDE_DIRS and
# HYPRE_LIBRARIES variables.
#
# The variables HYPRE_INSTALL_PREFIX or ENV HYPRE_ROOT control
# search locations for include and libraries.
#

set(hypre_search_paths ${HYPRE_DIR} ENV HYPRE_ROOT)
find_library(HYPRE_LIBRARY HYPRE HINTS ${hypre_search_paths} PATH_SUFFIXES lib)
find_path(HYPRE_INCLUDE_DIR HYPRE.h HINTS ${hypre_search_paths} PATH_SUFFIXES include)

if(HYPRE_INCLUDE_DIR)
  set(HYPRE_CONFIG_FILE ${HYPRE_INCLUDE_DIR}/HYPRE_config.h)
  include(SearchHeaderFile)
  search_header_file(${HYPRE_CONFIG_FILE} "HYPRE_HAVE_MPI" HYPRE_IS_PARALLEL)
  if(HYPRE_IS_PARALLEL)
    set(HYPRE_IS_PARALLEL True)
  else()
    set(HYPRE_IS_PARALLEL False)
  endif()
  if(NOT HYPRE_VERSION)
    search_header_file(${HYPRE_CONFIG_FILE} HYPRE_RELEASE_VERSION HYPRE_VERSION STRIP_QUOTES)
  endif()
  unset(HYPRE_CONFIG_FILE)
endif()

if(HYPRE_IS_PARALLEL)
  find_package(MPI)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE
				  REQUIRED_VARS HYPRE_LIBRARY HYPRE_INCLUDE_DIR
                                  VERSION_VAR HYPRE_VERSION)

if(HYPRE_FOUND)
  set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
  list(APPEND HYPRE_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
  if(HYPRE_IS_PARALLEL AND MPI_C_FOUND)
    list(APPEND HYPRE_LIBRARIES ${MPI_C_LIBRARIES})
    list(APPEND HYPRE_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH})
  endif()
  if(NOT TARGET hypre)
    add_library(hypre UNKNOWN IMPORTED)
    set_target_properties(hypre PROPERTIES
        IMPORTED_LOCATION "${HYPRE_LIBRARY}"
        #IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${HYPRE_LIBRARIES}")
        #INTERFACE_COMPILE_DEFINITIONS "${MPI_C_COMPILE_FLAGS}")
  endif()
endif()
