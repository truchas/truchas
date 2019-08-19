unset(CMAKE_Fortran_MODULE_DIRECTORY)


# Directory to install objects
# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.11...3.15)


# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(IRL VERSION 0.1
	LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)

set(IRL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/packages/irl")
set(IRL_BINARY_DIR "${PROJECT_BINARY_DIR}")
set(IRL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

if(CMAKE_BUILD_TYPE)
  if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    list(APPEND IRL_CXX_FLAGS "-DNDEBUG -DNDEBUG_PERF")
  endif()
endif()

# Inform CMake where to look for CMake files.
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_INSTALL_PREFIX}/cmake/Modules"
		      "${IRL_SOURCE_DIR}/cmake/Utils")

# Check valid build type
# Validate a user-supplied value for CMAKE_BUILD_TYPE.
if(CMAKE_BUILD_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)
  if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "^(DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be {Debug|Release|RelWithDebInfo|MinSizeRel}")
  endif()
endif()

# Directory to install objects
set(IRL_INSTALL_DIR "${IRL_BINARY_DIR}/install/${CMAKE_BUILD_TYPE}")
set(CMAKE_Fortran_MODULE_DIRECTORY "${IRL_INSTALL_DIR}/include")

# Find external packages needed by IRL
set(EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tpl/Eigen")
set(EIGEN_INCLUDE_DIR "${EIGEN_DIR}")


# Other common choices are STATIC, SHARED, and MODULE
# Output libname matches target name, with the usual extensions on your system
add_library(irl STATIC)
add_library(irl_c STATIC)
add_library(irl_fortran STATIC)

set_property(TARGET irl PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/libs")
set_property(TARGET irl PROPERTY CXX_STANDARD 14)
set_property(TARGET irl PROPERTY COMPILE_FLAGS "${IRL_CXX_FLAGS}")
set_property(TARGET irl_c PROPERTY CXX_STANDARD 14)
set_property(TARGET irl_c PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/libs")
set_property(TARGET irl_c PROPERTY COMPILE_FLAGS "${IRL_CXX_FLAGS}")

# C++ IRL Base
target_include_directories(irl PRIVATE ./)

# C Interface
target_include_directories(irl_c PRIVATE ./)
target_link_libraries(irl_c PRIVATE irl)

# Fortran Interface
set_target_properties(irl_fortran 
		      PROPERTIES FORTRAN_MODULE_DIRECTORY "${IRL_INSTALL_DIR}/include")
set_property(TARGET irl_fortran PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/libs")
set_property(TARGET irl_fortran PROPERTY COMPILE_FLAGS "${CMAKE_FORTRAN_FLAGS}")
target_include_directories(irl_fortran PRIVATE "${IRL_INSTALL_DIR}/include")
target_link_libraries(irl_fortran PRIVATE irl)
target_link_libraries(irl_fortran PRIVATE irl_c)

# Directory for test files
# Will automatically download GoogleTest and compile it to
# use for the testing framework.
#include(CTest)
#enable_testing()
#add_subdirectory(./tests)

# Need to find and link to Eigen package for C++ and C.
target_include_directories(irl SYSTEM PRIVATE "${EIGEN_INCLUDE_DIR}")
target_include_directories(irl_c SYSTEM PRIVATE "${EIGEN_INCLUDE_DIR}")

# Directory with source files.
add_subdirectory("${IRL_SOURCE_DIR}/src")

# Directory with development tools.
add_subdirectory("${IRL_SOURCE_DIR}/tools")

# Export targets to be used by other programs
install(TARGETS irl EXPORT irl_exp DESTINATION "${IRL_INSTALL_DIR}/libs")
install(TARGETS irl_c EXPORT irl_c_exp DESTINATION "${IRL_INSTALL_DIR}/libs")
install(TARGETS irl_fortran EXPORT irl_fortran_exp DESTINATION "${IRL_INSTALL_DIR}/libs")
install(EXPORT irl_exp DESTINATION "${IRL_INSTALL_DIR}/cmake") 
install(EXPORT irl_c_exp DESTINATION "${IRL_INSTALL_DIR}/cmake")
install(EXPORT irl_fortran_exp DESTINATION "${IRL_INSTALL_DIR}/cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY "${Truchas_MODULE_DIR}")
