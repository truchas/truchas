# This module finds the AmgX library. If successful, it defines the imported
# library target amgx. It is generally sufficient to include amgx as a target
# link library; CMake will automatically handle adding the needed compile
# include flags needed for the .mod files.
#
# To provide the module with hints on where to find your AMGX installations you
# have several options. You can include the root installation directories in the
# setting of the CMAKE_PREFIX_PATH variable, or you can set the environment
# variables or CMake variables AMGX_ROOT.
#
# This module uses a workaround for compiling with Darwin: it forces the gcc
# 5.5.0 stdc++ library into the compilation. This is because AmgX must be built
# with gcc, and it expects a standard library newer than the Darwin system's to
# be used.

find_library(AMGX_LIBRARY amgxsh
  HINTS ${AMGX_ROOT} ENV AMGX_ROOT
  PATH_SUFFIXES lib)
list(APPEND AMGX_LIBRARY "/projects/opt/centos7/gcc/5.5.0/lib64/libstdc++.so")
find_path(AMGX_INCLUDE_DIR amgx_c.h
  HINTS ${AMGX_ROOT} ENV AMGX_ROOT
  PATH_SUFFIXES include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMGX DEFAULT_MSG AMGX_INCLUDE_DIR AMGX_LIBRARY)

add_library(amgx INTERFACE IMPORTED)

set_target_properties(amgx PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${AMGX_INCLUDE_DIR}
  INTERFACE_LINK_LIBRARIES "${AMGX_LIBRARY}")
