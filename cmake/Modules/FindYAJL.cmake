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

# --- Function to handle the get_filename_component change in 2.8.11
function(get_directory var f)

  set(_comp_type)
  if ( ${CMAKE_VERSION} VERSION_LESS 2.8.11 )
    set(_comp_type PATH)
  else()
    set(_comp_type DIRECTORY)
  endif()
  get_filename_component(tmp ${f} ${_comp_type})
  set(${var} ${tmp} PARENT_SCOPE)

endfunction(get_directory var f)

# --- Search paths
set(_yajl_search_paths
    ${YAJL_INSTALL_PREFIX}
    ENV YAJL_ROOT)

# --- Locate include path
find_path(YAJL_INCLUDE_DIR
          NAMES yajl/yajl_common.h
          HINTS ${_yajl_search_paths}
          PATH_SUFFIXES include)

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
if ( YAJL_LIBRARY_STATIC OR YAJL_LIBRARY_SHARED )
  if ( YAJL_LIBRARY_STATIC )
    #get_filename_component(YAJL_LIBRARY_DIR ${YAJL_LIBRARY_STATIC} DIRECTORY)
    get_directory(YAJL_LIBRARY_DIR ${YAJL_LIBRARY_STATIC})
  endif()
  if ( YAJL_LIBRARY_SHARED )
    #get_filename_component(YAJL_LIBRARY_DIR ${YAJL_LIBRARY_SHARED} DIRECTORY)
    get_directory(YAJL_LIBRARY_DIR ${YAJL_LIBRARY_SHARED})
  endif()
else ()  
  set(YAJL_LIBRARY_DIR YAJL_LIBRARY_DIR-NOTFOUND)
endif()

# --- Define the version

# Do not repeat this search if already defined
if ( NOT DEFINED YAJL_VERSION )

  # Default to not found state
  set(YAJL_VERSION YAJL_VERSION-NOTFOUND)

  find_file(_yajl_version_h
            NAMES yajl/yajl_version.h
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

