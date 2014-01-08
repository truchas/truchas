#
# Find the yajl library and header files
#
# YAJL_FOUND          - true if yajl was found
#
# If found, the following variables are defined:
#
# YAJL_INCLUDE_DIR    - where to find the header file yajl/yajl_common.h
# YAJL_LIBRARY_DIR    - where to find the library (libyajl_s.a, libyajl.so)
# YAJL_LIBRARY        - the library to use
#
# YAJL_VERSION_MAJOR  - the major version of yajl
# YAJL_VERSION_MINOR  - the minor version of yajl
# YAJL_VERSION_PATCH  - the patch version of yajl
# YAJL_VERSION_STRING - the yajl version string (x.y.z)
#
# An includer may define YAJL_INCLUDE_DIR and YAJL_LIBRARY_DIR to tell this
# module where to look.
#

find_path(YAJL_INCLUDE_DIR yajl/yajl_version.h)

if(YAJL_INCLUDE_DIR AND EXISTS "${YAJL_INCLUDE_DIR}/yajl/yajl_version.h")
  file(STRINGS "${YAJL_INCLUDE_DIR}/yajl/yajl_version.h" YAJL_H REGEX "^#define YAJL_MAJOR .*$")
  string(REGEX REPLACE "^.*YAJL_MAJOR ([0-9]+).*$" "\\1" YAJL_VERSION_MAJOR "${YAJL_H}")
  file(STRINGS "${YAJL_INCLUDE_DIR}/yajl/yajl_version.h" YAJL_H REGEX "^#define YAJL_MINOR .*$")
  string(REGEX REPLACE "^.*YAJL_MINOR ([0-9]+).*$" "\\1" YAJL_VERSION_MINOR "${YAJL_H}")
  file(STRINGS "${YAJL_INCLUDE_DIR}/yajl/yajl_version.h" YAJL_H REGEX "^#define YAJL_MICRO .*$")
  string(REGEX REPLACE "^.*YAJL_MICRO ([0-9]+).*$" "\\1" YAJL_VERSION_PATCH "${YAJL_H}")
  set(YAJL_VERSION_STRING "${YAJL_VERSION_MAJOR}.${YAJL_VERSION_MINOR}.${YAJL_VERSION_PATCH}")
endif()
 
find_library(YAJL_LIBRARY yajl yajl_s HINTS ${YAJL_LIBRARY_DIR})
 
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAJL REQUIRED_VARS YAJL_LIBRARY YAJL_INCLUDE_DIR
                                  VERSION_VAR YAJL_VERSION_STRING)

if(YAJL_FOUND)
  get_filename_component(YAJL_LIBRARY_DIR ${YAJL_LIBRARY} DIRECTORY)
  set_target_properties(yajl PROPERTIES
                        IMPORTED_LOCATION ${YAJL_LIBRARY_DIR}
                        IMPORTED_INCLUDE_DIRECTORIES ${YAJL_INCLUDE_DIR})
endif()

