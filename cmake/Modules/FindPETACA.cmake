# - Find petaca
# ############################################################################ #
# Find the PETACA includes and library.
# Once done this will define
#
#  PETACA_MODULE_DIR     (PATH) - where to find the Fortran module files
#  PETACA_INCLUDE_DIRS   (PATH) - list of paths to include to compile
#  PETACA_LIBRARY        (FILE) - Petaca library
#  PETACA_LIBRARIES      (FILE) - list of libraries required to link against
#                                 Petaca
#  PETACA_FOUND          (BOOL) - True if Petaca found.
#
#  PETACA_VERSION        (STRING) - The version of Petaca found (x.y.z)
#
# PETACA will be added to the packages directory and this module file will be
# deprecated. For now, need this file to avoid rebuilding petaca when it is not
# required.
# Module file searches for the library file then builds the include directory
# assuming that /some/path/lib/libpetaca.a is the pattern for the library file
# and the include path will be /some/path/include. Search for petaca is
# controlled by PETACA_INSTALL_PREFIX, otherwise find_library looks in the
# standard CMake locations. See the cmake --help documentation for details.
#
# ############################################################################ #

if(PETACA_FIND_QUIETLY)
  set(_FIND_YAJL_ARG QUIET)
endif()
find_package(YAJL "2.0.4" ${_FIND_YAJL_ARG})

if(YAJL_FOUND)
  find_library(PETACA_LIBRARY petaca)

  # Module files are installed with the library file.
  if(PETACA_LIBRARY)
    get_filename_component(PETACA_MODULE_DIR ${PETACA_LIBRARY} DIRECTORY)
  endif()

  # No version number is currently available
  set(PETACA_VERSION PETACA_VERSION-NOTFOUND)

  if(PETACA_LIBRARY AND PETACA_MODULE_DIR)
    set(PETACA_LIBRARIES ${PETACA_LIBRARY} ${YAJL_LIBRARY})
    set(PETACA_INCLUDE_DIRS ${PETACA_MODULE_DIR})
    if(NOT TARGET petaca)
      add_library(petaca UNKNOWN IMPORTED)
      set_target_properties(petaca PROPERTIES
          IMPORTED_LOCATION "${PETACA_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${PETACA_INCLUDE_DIRS}"
          INTERFACE_LINK_LIBRARIES "${YAJL_LIBRARY}")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETACA REQUIRED_VARS PETACA_LIBRARY PETACA_MODULE_DIR)
