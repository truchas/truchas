# This module finds the Fortran FCUDA library. If successful, it defines the
# imported library target fcuda. It is generally sufficient to include fcuda
# as a target link library; CMake will automatically handle adding the compile
# include flags needed for the .mod files.
#
# To provide the module with hints on where to find your FCUDA
# installations you have several options. You can include the root installation
# directories in the setting of the CMAKE_PREFIX_PATH variable, or you can set
# the environment variables or CMake variables FCUDA_ROOT.

find_library(FCUDA_LIBRARY
  NAMES fcuda
  HINTS ${FCUDA_ROOT} ENV FCUDA_ROOT
  PATH_SUFFIXES lib)
find_path(FCUDA_INCLUDE_DIR
  NAMES fcuda.mod
  HINTS ${FCUDA_ROOT} ENV FCUDA_ROOT
  PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FCUDA DEFAULT_MSG FCUDA_INCLUDE_DIR FCUDA_LIBRARY)

add_library(fcuda INTERFACE IMPORTED)

set_target_properties(fcuda PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${FCUDA_INCLUDE_DIR}
  INTERFACE_LINK_LIBRARIES ${FCUDA_LIBRARY})
