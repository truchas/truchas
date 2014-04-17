# - Find yajl
# ############################################################################ #
# Find the YAJL includes and library.
# Once done this will define
#
#  YAJL_INCLUDE_DIR    (PATH) - where to find yajl/yajl_common.h, etc.
#  YAJL_LIBRARY_SHARED (FILE) - Shared yajl library
#  YAJL_LIBRARY_STATIC (FILE) - Static yajl library
#  YAJL_LIBRARY_DIR    (FILE) - Location of the yajl libraries
#  YAJL_FOUND          (BOOL) - True if yajl found.
#
#  YAJL_VERSION        (STRING) - The version of yajl found (x.y.z)

# ############################################################################ #

# --- Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# --- Search paths
set(_yajl_search_paths
    ${YAJL_INSTALL_PREFIX}
    ENV YAJL_ROOT)

# --- Locate include path
find_path(YAJL_INCLUDE_DIR
          NAMES yajl_common.h
          HINTS ${_yajl_search_paths}
          PATH_SUFFIXES include include/yajl)

# --- Locate the static library
find_library(YAJL_LIBRARY_STATIC
             NAMES yajl_s
             HINTS ${_yajl_search_paths}
             PATH_SUFFIXES lib)
        
# --- Locate the shared library        
find_library(YAJL_LIBRARY_SHARED
             NAMES yajl
             HINTS ${_yajl_search_paths}
             PATH_SUFFIXES lib)

# --- Define the library directory
find_path(YAJL_LIBRARY_DIR
          ${YAJL_LIBRARY_SHARED}
          HINTS ${_yajl_search_paths}
          PATH_SUFFIXES lib)

# --- Define the version

# Do not repeat this search if already defined
if ( NOT DEFINED YAJL_VERSION )

  # Default to not found state
  set(YAJL_VERSION YAJL_VERSION-NOTFOUND)

  find_file(_yajl_version_h
            NAMES yajl_version.h
            HINTS ${YAJL_INCLUDE_DIR})

  if ( _yajl_version_h )
    file(STRINGS "${_yajl_version_h}" _yajl_version_major_string 
                 REGEX "^#define YAJL_MAJOR")
    string(REGEX REPLACE "[^0-9]" "" _yajl_version_major "${_yajl_version_major_string}")            
  
    file(STRINGS "${_yajl_version_h}" _yajl_version_minor_string 
                 REGEX "^#define YAJL_MINOR")
    string(REGEX REPLACE "[^0-9]" "" _yajl_version_minor "${_yajl_version_minor_string}")            

    file(STRINGS "${_yajl_version_h}" _yajl_version_micro_string 
                 REGEX "^#define YAJL_MICRO")
    string(REGEX REPLACE "[^0-9]" "" _yajl_version_micro "${_yajl_version_micro_string}")            

    if ( DEFINED _yajl_version_major AND
         DEFINED _yajl_version_minor AND
         DEFINED _yajl_version_micro)
       set(YAJL_VERSION
         ${_yajl_version_major}.${_yajl_version_minor}.${_yajl_version_micro})
    endif()

  endif()

endif()  

# --- Print a useful status to the screen
find_package_handle_standard_args(YAJL
                                  REQUIRED_VARS YAJL_INCLUDE_DIR
                                  VERSION_VAR YAJL_VERSION
                                  HANDLE_COMPONENTS)

