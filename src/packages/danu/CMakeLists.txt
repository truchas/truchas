# ############################################################################ #
#
#   DANU  
#  
# ############################################################################ #


# ---------------------------------------------------------------------------- #
# Enabled interfaces (Default: OFF) 
# ---------------------------------------------------------------------------- #
set(ENABLE_Fortran ON)
set(ENABLE_Python  ON)

# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(Danu C CXX Fortran)

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH
    "${CMAKE_MODULE_PATH}"
    "${Danu_SOURCE_DIR}/cmake")
include(DanuGlobalMacros)

# project version
set( Danu_MAJOR_VERSION 1 )
set( Danu_MINOR_VERSION 0 )
set( Danu_PATCH_LEVEL 0 )
  set( Danu_VERSION "${Danu_MAJOR_VERSION}.${Danu_MINOR_VERSION}")
if ( ${Danu_PATCH_LEVEL} GREATER 0 )
  set( Danu_VERSION "${Danu_VERSION}-${Danu_PATCH_LEVEL}")
endif()
message(STATUS "Danu Version: ${Danu_VERSION}")


# ---------------------------------------------------------------------------- #
# Build Configuration
# ---------------------------------------------------------------------------- #

# Default build configuration is Release (optimization)
set(Danu_BUILD_TYPE_DFLT Release)
if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE ${Danu_BUILD_TYPE_DFLT})
endif()
message(STATUS "Build configuration type ${CMAKE_BUILD_TYPE}")

# Need the upper case name to report build flags at the end
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UC)

# Build the unit test suite 
option(ENABLE_Tests  "Build the Unit Test suite" ON)
if ( ENABLE_Tests )
  enable_testing()
endif()

# Build with PIC (Position-independent code)
option(ENABLE_PIC "Build position-independent code" OFF)
if (ENABLE_Python)
  set(ENABLE_PIC TRUE 
      CACHE BOOL "Force ENABLE_PIC to True since ENABLE_Python is TRUE" FORCE)
endif()

# ---------------------------------------------------------------------------- #
# Compiler Configuration 
# ---------------------------------------------------------------------------- #
include(CompilerConfig)

# ---------------------------------------------------------------------------- #
# Installation definitions 
# ---------------------------------------------------------------------------- #

# Binary files
if ( NOT Danu_BIN_INSTALL_DIR )
  set(Danu_BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
endif()

# Include files
if ( NOT Danu_INCLUDE_INSTALL_DIR)
  set(Danu_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Module (Fortran) files
if ( NOT Danu_MODULE_INSTALL_DIR)
  set(Danu_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Library files
if ( NOT Danu_LIB_INSTALL_DIR)
  set(Danu_LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Python modules
if ( ENABLE_Python )
  if ( NOT Danu_Python_INSTALL_DIR )
    set(Danu_Python_INSTALL_DIR 
        ${Danu_LIB_INSTALL_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
  endif()
endif()  

# CMake configuration files
if ( NOT Danu_CONFIG_INSTALL_DIR)
  set(Danu_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/danu)
endif()

# ---------------------------------------------------------------------------- #
# External Software Package Configuration
# ---------------------------------------------------------------------------- #
#DEBUG
#DEBUG# HDF5 (Required)
#DEBUGfind_package(HDF5)
#DEBUGif ( NOT HDF5_FOUND )
#DEBUG
#DEBUG  # Set the HDF5 build flag
#DEBUG  set(BUILD_HDF5 TRUE)
#DEBUG  message(STATUS "Will build HDF5")
#DEBUG
#DEBUG  # Look for zlib too
#DEBUG  find_package(ZLIB)
#DEBUG  if (NOT ZLIB_FOUND OR ( ${ZLIB_VERSION_STRING} VERSION_LESS "1.2.8" ) )
#DEBUG    set(BUILD_ZLIB TRUE)
#DEBUG    message(STATUS "Could not locate a sufficient zlib will build")
#DEBUG  endif()
#DEBUG
#DEBUG
#DEBUGendif()    
#DEBUG
#DEBUG
#DEBUG# Search for Python, NumPy and SWIG
#DEBUGif ( ENABLE_Python )
#DEBUG
#DEBUG  # Python (REQUIRED)
#DEBUG  find_package(Python)
#DEBUG  if ( NOT PYTHON_FOUND )
#DEBUG    message(FATAL_ERROR "Enabled Python interface but could not locate "
#DEBUG                        "Python. Specify a python executable with "
#DEBUG			"-D PYTHON_EXECUTABLE:STRING=")
#DEBUG  endif()
#DEBUG
#DEBUG  # NumPy (REQUIRED)
#DEBUG  find_package(NumPy)
#DEBUG  if ( NOT NUMPY_FOUND )
#DEBUG    message(FATAL_ERROR "The ${Project} Python interface requires the "
#DEBUG                        "Python package NumPy (Numerical Python). "
#DEBUG			"Please install this package and re-run cmake.")
#DEBUG  endif()
#DEBUG
#DEBUG  # SWIG
#DEBUG  find_package(SWIG)
#DEBUG  if ( NOT SWIG_FOUND OR ( "${SWIG_VERSION}" VERSION_LESS "2.0" ) )
#DEBUG    message(STATUS "Can not locate sufficient SWIG package. Add SWIG to build")
#DEBUG    set(BUILD_SWIG TRUE)
#DEBUG  endif()  
#DEBUG
#DEBUGendif()
#DEBUG
#DEBUG# Unit tests
#DEBUGif ( ENABLE_Tests )
#DEBUG
#DEBUG  find_package(Check)
#DEBUG
#DEBUG  if ( NOT CHECK_FOUND )
#DEBUG    set(BUILD_Check TRUE)
#DEBUG    message(STATUS "Will build Check for the test suite")
#DEBUG  endif()
#DEBUG
#DEBUGendif()  
#DEBUG
#DEBUG
#DEBUG# Build the external packages
#DEBUGinclude(BuildExternalPackages)
#DEBUG
# ---------------------------------------------------------------------------- #
# Danu Source Directory 
# ---------------------------------------------------------------------------- #
add_subdirectory(src)

#DEBUG# ---------------------------------------------------------------------------- #
#DEBUG# Testing 
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUGif(ENABLE_Tests)
#DEBUG  message(STATUS "Build the unit test suite")
#DEBUG  add_subdirectory(test)
#DEBUGendif()  
        
# ---------------------------------------------------------------------------- #
# External configuration files 
# ---------------------------------------------------------------------------- #
# Create the target file
set(Danu_TARGETS_FILE ${Danu_BINARY_DIR}/danu-targets.cmake)
set(Danu_TARGETS danu)
if(ENABLE_Fortran)
  list(APPEND Danu_TARGETS danu_f)
endif()
export(TARGETS ${Danu_TARGETS} FILE ${Danu_TARGETS_FILE})


# Create danu-config.cmake (find_package(Danu))
set(Danu_CONFIG_FILE ${Danu_BINARY_DIR}/danu-config.cmake)
configure_file(${Danu_SOURCE_DIR}/danu-config.in
               ${Danu_CONFIG_FILE}
	       @ONLY)
install(FILES ${Danu_CONFIG_FILE} ${Danu_TARGETS_FILE} 
        DESTINATION ${Danu_CONFIG_INSTALL_DIR})



#DEBUG# ---------------------------------------------------------------------------- #
#DEBUG# Print out final configuration
#DEBUG# ---------------------------------------------------------------------------- #
#DEBUGmessage(STATUS) 
#DEBUGmessage(STATUS "-------------------------------------------------------") 
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "Installation prefix ${CMAKE_INSTALL_PREFIX}")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "Build Type      = ${CMAKE_BUILD_TYPE}")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "C COMPILER      = ${CMAKE_C_COMPILER}")
#DEBUGmessage(STATUS "CMAKE_C_FLAG\t=${CMAKE_C_FLAGS}")
#DEBUGmessage(STATUS "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC}\t"
#DEBUG               "=${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "Fortran enabled = ${ENABLE_Fortran}")
#DEBUGif(ENABLE_Fortran)
#DEBUG    message(STATUS "Fortran COMPILER=${CMAKE_Fortran_COMPILER}")
#DEBUG    message(STATUS "CMAKE_Fortran_FLAG\t=${CMAKE_Fortran_FLAGS}")
#DEBUG    message(STATUS "CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}\t"
#DEBUG                     "=${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
#DEBUGendif()
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "Python enabled  = ${ENABLE_Python}")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "External software packages")
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "HDF5 Found = ${HDF5_FOUND}")
#DEBUGmessage(STATUS "HDF5 INCLUDE DIRS = ${HDF5_INCLUDE_DIRS}")
#DEBUGmessage(STATUS "HDF5 LIBRARIES = ${HDF5_LIBRARIES}")
#DEBUGmessage(STATUS)
#DEBUGif(ENABLE_Tests)
#DEBUG  message(STATUS)
#DEBUG  message(STATUS "Check Found = ${Check_FOUND}")
#DEBUG  message(STATUS "Check INCLUDE DIRS = ${Check_INCLUDE_DIRS}")
#DEBUG  message(STATUS "Check LIBRARIES = ${Check_LIBRARIES}")
#DEBUG  message(STATUS)
#DEBUGendif()  
#DEBUGif(ENABLE_Python)
#DEBUG  message(STATUS "")
#DEBUG  message(STATUS "SWIG=${SWIG_EXECUTABLE}")
#DEBUG  message(STATUS "SWIG VERSION=${SWIG_VERSION}")
#DEBUG  message(STATUS "")
#DEBUG  message(STATUS "Python Executable=${PYTHON_EXECUTABLE}")
#DEBUG  message(STATUS "Python Version=${PYTHON_VERSION}")
#DEBUG  message(STATUS "Python Libraries=${PYTHON_LIBRARIES}")
#DEBUG  message(STATUS "Python Include Directories=${PYTHON_INCLUDE_DIRS}")
#DEBUG  message(STATUS "NumPy Version=${NUMPY_VERSION}")
#DEBUG  message(STATUS "NumPy Include Directories=${NUMPY_INCLUDE_DIRS}")
#DEBUG  message(STATUS "")
#DEBUGendif(ENABLE_Python)
#DEBUGmessage(STATUS)
#DEBUGmessage(STATUS "-------------------------------------------------------") 
#DEBUGmessage(STATUS)
