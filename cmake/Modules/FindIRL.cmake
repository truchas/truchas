# - Find IRL (Interface Reconstruction Library)
# ############################################################################ #
# Find the IRL includes and library.
# Once done this will define
#
#  IRL_INCLUDE_DIRS   (PATH) - list of paths to include to compile
#  IRL_LIBRARY        (FILE) - IRL library needed for use in Fortran
#  IRL_LIBRARY	     (FILES) - All IRL libraries required for linking
#
#  IRL_FOUND          (BOOL) - True if IRL found.
#
#  IRL_VERSION      (STRING) - The version of IRL found (x.y.z)
#
# Module file searches for the library file then builds the include directory
# assuming that /some/path/lib/libirl.a is the pattern for the library file
# and the include path will be /some/path/include. Search for IRL is
# controlled by IRL_INSTALL_PREFIX, otherwise find_library looks in the
# standard CMake locations. See the cmake --help documentation for details.
#
# ############################################################################ #

find_library(IRL_LIBRARY 
	     NAMES irl_fortran 
	     PATHS ENV IRL_LIB_PATH
	     NO_DEFAULT_PATH)

if(IRL_LIBRARY)

  # No version number is currently available
  set(IRL_VERSION IRL_VERSION-NOTFOUND)

  get_filename_component(IRL_LIB_DIR ${IRL_LIBRARY} DIRECTORY)
  set(IRL_DIR "${IRL_LIB_DIR}/..")
  set(BUILD_TESTING "OFF")
  add_subdirectory("${IRL_DIR}/abseil-cpp/" "${Truchas_BINARY_DIR}/abseil-cpp")
  set(IRL_LIBRARIES ${IRL_LIB_DIR}/libirl.a ${IRL_LIB_DIR}/libirl_c.a ${IRL_LIBRARY} absl::flat_hash_map absl::inlined_vector -lstdc++)
  set(IRL_INCLUDE_DIRS "${IRL_DIR}/include/")
  if(NOT TARGET irl)
    add_library(irl UNKNOWN IMPORTED)
    set_target_properties(irl PROPERTIES
        IMPORTED_LOCATION "${IRL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${IRL_INCLUDE_DIRS}"
	INTERFACE_LINK_LIBRARIES "${IRL_LIBRARIES}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IRL 
				  REQUIRED_VARS IRL_LIBRARIES IRL_INCLUDE_DIRS
				  VERSION_VAR IRL_VERSION)
