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

# --- Locate Petaca library
find_library(PETACA_LIBRARY
              NAMES petaca
              HINTS ${PETACA_INSTALL_PREFIX}
              PATH_SUFFIXES lib)

# --- Define the include path from the library file
if ( PETACA_LIBRARY ) 
  get_directory(_petaca_lib_dir ${PETACA_LIBRARY})
  get_directory(PETACA_MODULE_DIR ${_petaca_lib_dir})
  set(PETACA_MODULE_DIR ${PETACA_MODULE_DIR}/include)

  if ( NOT EXISTS ${PETACA_MODULE_DIR} ) 
    message(WARNING "Failed to define PETACA_MODULE_DIR")
    set(PETACA_MODULE_DIR PETACA_MODULE_DIR-NOTFOUND)
  endif()
  
endif()

# --- Define the version, hard coded for now
if ( PETACA_LIBRARY AND PETACA_INCLUDE_DIR )
  set(PETACA_VERSION 2ba5af4)
else()  
  set(PETACA_VERSION PETACA_VERSION-NOTFOUND)
endif()

# --- Search for YAJL 
find_package(YAJL)
if ( NOT YAJL_FOUND ) 
  message(WARNING "Could not locate YAJL, required package for PETACA")
endif()

# --- Define PETACA_INCLUDE_DIRS
set(PETACA_INCLUDE_DIRS ${PETACA_MODULE_DIR})
if ( YAJL_FOUND AND PETACA_INCLUDE_DIRS)
  list(APPEND PETACA_INCLUDE_DIRS ${YAJL_INCLUDE_DIR})
  list(REMOVE_DUPLICATES PETACA_INCLUDE_DIRS)
endif()

# --- Define PETACA_LIBRARIES
set(PETACA_LIBRARIES ${PETACA_LIBRARY})
if ( YAJL_FOUND  AND PETACA_LIBRARY)
  if ( BUILD_SHARED_LIBS )
    list(APPEND PETACA_LIBRARIES ${YAJL_LIBRARY_SHARED})
  else()
    list(APPEND PETACA_LIBRARIES ${YAJL_LIBRARY_STATIC})
  endif()  
endif()


# --- Print a useful status to the screen
find_package_handle_standard_args(PETACA
                                  REQUIRED_VARS 
                                    PETACA_LIBRARY
                                    PETACA_MODULE_DIR
                                  VERSION_VAR 
                                    PETACA_VERSION
                                  HANDLE_COMPONENTS)

