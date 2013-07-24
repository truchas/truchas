# ############################################################################ #
# 
# SPLIT_LIBRARY_NAME(<file>|target
#                    PATH path_var
#                    FILE filename_var
#                    LIBNAME libname_var
#                    EXT ext_var)
#
#
# Given the full path name for a library or an import library target
# split the library file name into PATH, LIBNAME and extension.
#
# Example:
#   SPLIT_LIBRARY_NAME(/usr/local/lib/librt.so
#                     PATH rt_path
#                     LIBNAME rt_libname
#                     EXT rt_ext)
# The variables rt_path, rt_libname and rt_ext would be
# set to:
#  rt_path=/usr/local/lib
#  rt_libname=rt
#  rt_ext=so
#
# This information can then used to build a -L/-l flag options
# -L${rt_path} -lrt
# on UNIX platforms.
include(CMakeParseArguments)
MACRO(SPLIT_LIBRARY_NAME library)

  set(options "")
  set(oneValue "PATH;LIBNAME;EXT")
  set(multiValue "")
  cmake_parse_arguments(PARSE "${options}" "${oneValue}" "${multiValue}" ${ARGN})

  if (TARGET "${library}")
    include(CMakeExpandImportedTargets)
    cmake_expand_imported_targets(_fullname LIBRARIES ${library})
  else()
    get_filename_component(_fullname "${library}" ABSOLUTE)
  endif()

  get_filename_component(${PARSE_PATH} "${_fullname}" PATH)
  get_filename_component(_filename "${_fullname}" NAME_WE)
  get_filename_component(${PARSE_EXT} "${_fullname}" EXT)

  set(lib_prefix_pattern)
  if ( "${${PARSE_EXT}}" MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(lib_prefix_pattern ${CMAKE_STATIC_LIBRARY_PREFIX})
  elseif("${${PARSE_EXT}}" MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(lib_prefix_pattern ${CMAKE_SHARED_LIBRARY_PREFIX})
  else()
    message(FATAL_ERROR "Unknown library suffix ${${PARSE_EXT}}")
  endif()

  string(REPLACE "${lib_prefix_pattern}" "" ${PARSE_LIBNAME} "${_filename}")

  #print_variable(PARSE_PATH)
  #print_variable(${PARSE_PATH})
  #print_variable(PARSE_LIBNAME)
  #print_variable(${PARSE_LIBNAME})
  #print_variable(PARSE_EXT)
  #print_variable(${PARSE_EXT})


ENDMACRO(SPLIT_LIBRARY_NAME library)
