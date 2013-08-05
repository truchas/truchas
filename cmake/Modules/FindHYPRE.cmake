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
include(FindPackageHandleStandardArgs)
include(PrintVariable)

# --- MACROS/FUNCTIONS
FUNCTION(_search_hypre_config_file config_h defmacro outvar)
  if ( EXISTS "${config_h}" )
    set(_regexp ".define[ \\t]+${defmacro}")
    set(_quote "\"")
    FILE(STRINGS "${config_h}" _tmp_string REGEX "${_regexp}")
    STRING(REGEX REPLACE ".define ${defmacro}[ \\t]+(.*)$" "\\1" _result "${_tmp_string}")
    STRING(REPLACE "${_quote}" " " _result "${_result}")
    STRING(REGEX REPLACE "[ \t\n]+" "" _result "${_result}")
  else()
    set(_result ${outvar}-NOTFOUND)
  endif()
  set(${outvar} ${_result} PARENT_SCOPE)

ENDFUNCTION()

# --- Begin main

# Search for HYPRE.h
find_path(HYPRE_INCLUDE_DIR
          NAMES HYPRE.h
          HINTS ENV HYPRE_ROOT ${HYPRE_INSTALL_PREFIX}
	  PATH_SUFFIXES include)

# Search for HYPRE library
find_library(HYPRE_LIBRARY
             NAMES HYPRE
             HINTS ENV HYPRE_ROOT ${HYPRE_INSTALL_PREFIX}
	     PATH_SUFFIXES lib)

# Define the cofig header file, will mine this to determine version
# and build type (MPI, non-MPI)
if ( HYPRE_INCLUDE_DIR )
  set(HYPRE_CONFIG_FILE ${HYPRE_INCLUDE_DIR}/HYPRE_config.h)
endif()

# Define the version
if ( NOT HYPRE_VERSION  AND HYPRE_CONFIG_FILE )
  _search_hypre_config_file(${HYPRE_CONFIG_FILE} HYPRE_RELEASE_VERSION HYPRE_VERSION)
endif()

# Define the build (MPI or non-MPI)
if (HYPRE_CONFIG_FILE)
  _search_hypre_config_file(${HYPRE_CONFIG_FILE} HYPRE_HAVE_MPI HYPRE_IS_PARALLEL)
  if(HYPRE_IS_PARALLEL)
    set(HYPRE_IS_PARALLEL True)
  else()
    set(HYPRE_IS_PARALLEL False)
  endif()
  if (HYPRE_IS_PARALLEL)
    find_package(MPI)
  endif()  
endif()

# Define the include directories, append MPI if needed
if ( NOT HYPRE_INCLUDE_DIRS)
  set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
  if (HYPRE_IS_PARALLEL)
    if (MPI_C_FOUND)
      list(APPEND HYPRE_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH})
    endif()
  endif()
endif()

# Define the  libraries
if ( NOT HYPRE_LIBRARIES )
  set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
  if (HYPRE_IS_PARALLEL)
    list(APPEND HYPRE_LIBRARIES ${MPI_C_LIBRARIES})
  endif()
endif()  

# Send a useful message if everything is found
find_package_handle_standard_args(HYPRE 
                                  VERSION_VAR HYPRE_VERSION
				  REQUIRED_VARS HYPRE_LIBRARY HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)

mark_as_advanced(
                 HYPRE_INCLUDE_DIR
                 HYPRE_CONFIG_FILE
                 HYPRE_LIBRARY)

