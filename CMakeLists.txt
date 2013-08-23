# ############################################################################ #
#                                                                              #
# Truchas                                                                      #
#                                                                              #
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Set policy and CMake version requirements
# ---------------------------------------------------------------------------- #
cmake_minimum_required( VERSION 2.8.8 )

# allow more human readable "if then else" constructs
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

if("${CMAKE_PATCH_VERSION}" GREATER 3 )
  cmake_policy(SET CMP0017 NEW)
endif()

# ---------------------------------------------------------------------------- #
# Project definitions
# ---------------------------------------------------------------------------- #
project(Truchas C CXX Fortran)

# Version
set(Truchas_MAJOR_VERSION 2)
set(Truchas_MINOR_VERSION 8)
set(Truchas_PATCH_VERSION dev)
set(Truchas_VERSION ${Truchas_MAJOR_VERSION}.${Truchas_MINOR_VERSION}.${Truchas_PATCH_VERSION})

# Copyright file
set(Truchas_COPYRIGHT_FILE ${Truchas_SOURCE_DIR}/COPYRIGHT)

# License file
set(Truchas_LICENSE_FILE ${Truchas_SOURCE_DIR}/LICENSE)

# Build date
try_run(BUILD_DATE_FAILED BUILD_DATE_COMPILED
        ${Truchas_BINARY_DIR} ${Truchas_SOURCE_DIR}/cmake/build_date.c
        RUN_OUTPUT_VARIABLE Truchas_BUILD_DATE)

if (BUILD_DATE_FAILED)
  message(WARNING "Failed to determine the build date stamp")
  if (NOT BUILD_DATE_COMPILED)
    message(WARNING "Could not compile build_date.c")
  endif()
  set(Truchas_BUILD_DATE NOW)
endif()

# Build hostname
try_run(HOSTNAME_FAILED HOSTNAME_COMPILED
        ${Truchas_BINARY_DIR} ${Truchas_SOURCE_DIR}/cmake/get_hostname.c
        RUN_OUTPUT_VARIABLE Truchas_HOSTNAME)

if (HOSTNAME_FAILED)
  message(WARNING "Failed to define hostname")
  if(NOT HOSTNAME_COMPILED)
    message(WARNING "Could not compile get_hostname.c")
  endif()
  set(Truchas_HOSTNAME ${CMAKE_HOST_SYSTEM_NAME}-generic)
endif()

# Location of project defined include and FindXXX.cmake files
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" 
                      "${Truchas_SOURCE_DIR}/cmake"
                      "${Truchas_SOURCE_DIR}/cmake/Modules"
                      "${Truchas_SOURCE_DIR}/cmake/Utils")

# Create all the Fortran module files in a directory at the root build directory
set(Truchas_MODULE_DIR ${Truchas_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${Truchas_MODULE_DIR})

# Modules/Compiler/NAG-Fortran.cmake is missing this.
if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set(CMAKE_Fortran_COMPILE_OPTIONS_PIC "-PIC")
endif()


# ---------------------------------------------------------------------------- #
# Global CMake Utilities 
# ---------------------------------------------------------------------------- #

# Useful print macro
include(PrintVariable)

# Set macros global_set, etc.
include(SetMacros)

# ---------------------------------------------------------------------------- #
# Build Options 
# ---------------------------------------------------------------------------- #

# Build type default:Release (optimized)
if ( NOT CMAKE_BUILD_TYPE )
  message(STATUS "CMAKE_BUILD_TYPE not set. Will use default Release (optimized)")
  set(CMAKE_BUILD_TYPE Release)
endif()    

# Build with MPI
option(ENABLE_MPI "Build the parallel (MPI) library" True)

# Build with PGSLib
option(ENABLE_PGSLib "Build with PGSLib support"  True)

# Build with UbikSolve
option(ENABLE_UbikSolve "Build with UbikSolve support"  True)

# Build with Tbrook 
option(ENABLE_TBROOK "Build with TBROOK support"  True)

# Build with Danu
option(ENABLE_Danu "Build with Danu support"  True)

# Build with Chaco
option(ENABLE_Chaco "Build with Chaco support"  True)

# Build with Chaparral
option(ENABLE_Chaparral "Build with Chaparral support"  True)

# Build with PIC (position in code)
option(WITH_PIC "Build with PIC support" True)

# Build shared libs
option(ENABLE_SHARED "Build shared libraries" False)

# Enable search for external packages
option(ENABLE_ExternalSearch "Enable external (NetCDF,HDF5,etc.) search" True)

# Test - need a swith here to enable testing
enable_testing()

# Set the BUILD_SHARED_LIBS variable
if(ENABLE_SHARED)
  set(BUILD_SHARED_LIBS True)
endif()

# Define the build id
include(GetTruchasBuildId)
get_truchas_build_id(Truchas_BUILD_ID)

# ---------------------------------------------------------------------------- #
# Fortran/C Interface 
# ---------------------------------------------------------------------------- #
include(FortranCInterface)

# Create the FortranCInterface Header (name mangling)
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface_names.h 
                         MACRO_NAMESPACE TR_ROUTINE_)
set(Truchas_FCIface_INCLUDE_DIR ${Truchas_BINARY_DIR})		       

# Check the compatibility of the Fortran/C compilers
# Test does fail when it shouldn't, hence QUIET
FortranCInterface_VERIFY(QUIET)

# Handle errors here
if ( NOT FortranCInterface_GLOBAL_FOUND )
  message(FATAL_ERROR "Could not find the global name-mangling for Fortran/C")
endif(NOT FortranCInterface_GLOBAL_FOUND )


# ---------------------------------------------------------------------------- #
# Setting Some Useful Language Variables 
# ---------------------------------------------------------------------------- #

# Because I get tried typing CMAKE_<LANG>_COMPILER_ID MATHCHES blah
# Set all to FALSE initially
foreach ( id GNU G95 INTEL NAG PGI ABSOFT )
  set(var_name Fortran_COMPILER_IS_${id})
  set(${var_name} FALSE)
endforeach()

# Now set the correct variable to true
string(TOUPPER "${CMAKE_Fortran_COMPILER_ID}" uc_id)
set(Fortran_COMPILER_IS_${uc_id} TRUE)
# ---------------------------------------------------------------------------- #
# Configuration Tests 
# ---------------------------------------------------------------------------- #

# newunit
try_compile(SUPPORTS_NEWUNIT
            ${Truchas_BINARY_DIR}
            ${Truchas_SOURCE_DIR}/cmake/newunit.f90)

# PIC support
if ( WITH_PIC OR ENABLE_TBROOK )
  include(CheckCCompilerFlag)
  include(UpdateCompilerFlags)
  set(test_pic_flags -fPIC -fpic -PIC -pic)
  set(pic_flag pic_flag-NOTFOUND)
  foreach(flag ${test_pic_flags})
    set(res_var pic_flag_ctest${flag})
    check_c_compiler_flag(${flag} ${res_var})
    if ( ${${res_var}} )
      set(pic_flag ${flag})
      break()
    endif(${${res_var}})
  endforeach()

  if(pic_flag)
    message(STATUS "Adding ${pic_flag} to CMAKE_C_FLAGS")
    update_compiler_flags(C FLAGS ${pic_flag})
  else()
    message(FATAL_ERROR "Could not determine the PIC flag for this C compiler. "
                        "Tried: ${test_pic_flags} all failed")
  endif()

  include(CheckCXXCompilerFlag)
  set(test_pic_flags -fPIC -fpic -PIC -pic)
  set(pic_flag pic_flag-NOTFOUND)
  foreach(flag ${test_pic_flags})
    set(res_var pic_flag_cxxtest${flag})
    check_cxx_compiler_flag(${flag} ${res_var})
    if ( ${${res_var}} )
      set(pic_flag ${flag})
      break()
    endif(${${res_var}})
  endforeach()

  if(pic_flag)
    message(STATUS "Adding ${pic_flag} to CMAKE_CXX_FLAGS")
    update_compiler_flags(CXX FLAGS ${pic_flag})
  else()
    message(FATAL_ERROR "Could not determine the PIC flag for this C++ compiler. "
                        "Tried: ${test_pic_flags} all failed")
  endif()

  # CMake does not have equivalent CheckFortran*.cmake modules
  # for now we'll hard code these flags by compiler type.
  set(pic_flag pic_flag-NOTFOUND)
  if(CMAKE_Fortran_COMPILE_OPTIONS_PIC)
    set(pic_flag ${CMAKE_Fortran_COMPILE_OPTIONS_PIC})
  else()
    if (Fortran_COMPILER_IS_NAG)
      set(pic_flag -PIC)
    elseif(Fortran_COMPILER_IS_INTEL)
      set(pic_flag -fpic)
    else()
      message(STATUS "PIC (position in code) compiler flag is not set for ${CMAKE_Fortran_COMPILER+ID}")
    endif()
  endif()

  if(pic_flag)
    message(STATUS "Adding ${pic_flag} to CMAKE_Fortran_FLAGS")
    update_compiler_flags(Fortran FLAGS ${pic_flag})
  else()
    message(FATAL_ERROR "Could not determine the PIC flag for this Fortran compiler.")
  endif()
  
endif()


# ---------------------------------------------------------------------------- #
# External Dependencies 
# ---------------------------------------------------------------------------- #

# MPI
if (ENABLE_MPI)

  # Attempt to find MPI without hints
  find_package(MPI)

  if ( NOT MPI_C_FOUND  OR NOT MPI_Fortran_FOUND)

    # Notify the user that the default search failed 
    message(STATUS "Initial MPI search failed")
    if (MPI_ROOT)
      message(STATUS "Search for MPI installation in ${MPI_ROOT}")
    endif()

    # First search failed, now look for mpicc in the usual places
    find_program(MPI_C_COMPILER
                 name mpicc
		 HINTS ${MPI_ROOT} ENV MPI_ROOT ENV MPIROOT ENV MPICH_ROOT
		 PATH_SUFFIXES bin)

    # First search failed, now look for mpicc in the usual places
    find_program(MPI_Fortran_COMPILER
                 NAMES mpif90 mpif77
		 HINTS ${MPI_ROOT} ENV MPI_ROOT ENV MPIROOT ENV MPICH_ROOT
		 PATH_SUFFIXES bin)

    # Attempt to locate the mpiexec command if the test suite is active
    if (ENABLE_TESTS)
      find_program(MPIEXEC
                   NAMES mpiexec mpirun aprun poe
	   	   HINTS ${MPI_ROOT} ENV MPI_ROOT ENV MPIROOT ENV MPICH_ROOT
		   PATH_SUFFIXES bin)
      if (NOT MPIEXEC_NUMPROC_FLAG)
	set(MPIEXEC_NUMPROC_FLAG -n)
      endif()

    endif()

    if (MPI_C_COMPILER OR MPI_Fortran_COMPILER)
      message(STATUS "Setting MPI_C_COMPILER=${MPI_C_COMPILER} "
	             " and MPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} to find MPI")
      find_package(MPI)
      if (NOT MPI_C_FOUND)
	message(WARNING "After setting MPI_C_COMPILER still can not locate MPI")
      endif()
    else()
      message(WARNING "Can not locate a C wrapper for MPI could not repeat MPI search")
    endif()

    # And now die because we really can not find MPI
    if (NOT MPI_C_FOUND)
      message(FATAL_ERROR "Can not locate an MPI installation"
                          "Either define an MPI installation with -D MPI_ROOT or"
                          "Explicitly set -D MPI_C_INCLUDE_PATH and "
			  "-D MPI_C_LIBRARIES to bypass this search")
    endif()			

  endif()  

endif()

# Python
find_package(Python)

# NumPy - Numerical Python
find_package(NumPy)

# LAPACK
find_package(LAPACK)

# BLAS
find_package(BLAS)




# ---------------------------------------------------------------------------- #
# Preprocessor definitions 
# ---------------------------------------------------------------------------- #

# Preprocessor - allow this to be a list to define global flags easily
if ( NOT Truchas_PREPROCESSOR )
  set(Truchas_PREPROCESSOR cpp -ansi)
endif()

# Turn off assertions in Truchas Fortran code
if ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" )
  list(APPEND Truchas_PREPROCESSOR -DNDEBUG)
endif()

# Add the compiler type to the global defines and the preprocessor
string(TOUPPER "-D${CMAKE_Fortran_COMPILER_ID}_COMPILER" Truchas_Compiler_Define)
if ( Fortran_COMPILER_IS_NAG )
  add_definitions(${Truchas_Compiler_Define})
  list(APPEND Truchas_PREPROCESSOR 
             ${Truchas_Compiler_Define})
else()
  add_definitions(${Truchas_Compiler_Define} ${Truchas_Compiler_Define}_WORKAROUND)
  list(APPEND Truchas_PREPROCESSOR 
              ${Truchas_Compiler_Define}
              ${Truchas_Compiler_Define}_WORKAROUND)
endif()

# Add USE_PGSLIB
if (ENABLE_PGSLib)
  list(APPEND Truchas_PREPROCESSOR -DUSE_PGSLIB)
endif()

# UbikSolve defines
if (ENABLE_UbikSolve)
  list(APPEND Truchas_PREPROCESSOR -DUSE_UBIKSOLVE)
endif()

# Chaco defines
if (ENABLE_Chaco)
  list(APPEND Truchas_PREPROCESSOR -DUSE_CHACO)
endif()

# Chaparral defines
if (ENABLE_Chaparral)
  list(APPEND Truchas_PREPROCESSOR -DCHAPARRAL_SUPPORT)
endif()

# NEWUNIT support
if (SUPPORTS_NEWUNIT)
  list(APPEND Truchas_PREPROCESSOR -DSUPPORTS_NEWUNIT)
endif()

# Arch x86_64 flag hard code for now
list(APPEND Truchas_PREPROCESSOR -DM64)

# OS defines -- hard code for now
list(APPEND Truchas_PREPROCESSOR -DLINUX)

# TBrook defines
if (ENABLE_TBROOK)
  list(APPEND Truchas_PREPROCESSOR -DUSE_TBROOK)
endif()

# Danu defines
if (ENABLE_Danu)
  list(APPEND Truchas_PREPROCESSOR -DUSE_DANU)
endif()

# More hardcoded defines - only used in drivers.f90 need to move
list(APPEND Truchas_PREPROCESSOR -DCODE_NAME='DUMMY')
list(APPEND Truchas_PREPROCESSOR -DUBIKSOLVE='UbikSolve')
list(APPEND Truchas_PREPROCESSOR -DPGSLIB='PGSLib')
list(APPEND Truchas_PREPROCESSOR -DCHACO='Chaco')
list(APPEND Truchas_PREPROCESSOR -DVERSION='${Truchas_VERSION}')
list(APPEND Truchas_PREPROCESSOR -DBUILD_DATE='${Truchas_BUILD_DATE}')
list(APPEND Truchas_PREPROCESSOR -DHOST_NAME='${Truchas_HOSTNAME}')
list(APPEND Truchas_PREPROCESSOR -DARCHITECTURE='${CMAKE_SYSTEM_PROCESSOR}')
list(APPEND Truchas_PREPROCESSOR -DFFLAGS_TRIMMED_1='${CMAKE_Fortran_FLAGS}')
list(APPEND Truchas_PREPROCESSOR -DFFLAGS_TRIMMED_2='${CMAKE_Fortran_FLAGS_RELEASE}')

# Add any user defined preprocessor flags
if (Truchas_FPP_FLAGS)
  list(APPEND Truchas_PREPROCESSOR ${Truchas_FPP_FLAGS})
endif()



# ---------------------------------------------------------------------------- #
# Installation Definitions 
# ---------------------------------------------------------------------------- #

# Need the paths defined before adding src to the build

# Default installation location is the Truchas_SOURCE_DIR
set(Truchas_DFLT_INSTALL_PREFIX ${Truchas_SOURCE_DIR}/install/${Truchas_BUILD_ID})
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${Truchas_DFLT_INSTALL_PREFIX}
      CACHE PATH "Installation prefix" FORCE)
endif()

# Binary installation path, override with -D Truchas_BIN_INSTALL_DIR
if (NOT Truchas_BIN_INSTALL_DIR)
  set(Truchas_BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
endif()

# Library installation path, override with -D Truchas_LIBRARY_INSTALL_DIR
if (NOT Truchas_LIBRARY_INSTALL_DIR)
  set(Truchas_LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Module installation path, override with -D Truchas_MODULE_INSTALL_DIR
if (NOT Truchas_MODULE_INSTALL_DIR)
  set(Truchas_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

# Copy entire include directory to  Truchas_MODULE_INSTALL_DIR
install(CODE "MESSAGE(STATUS \"Installing: modules in ${Truchas_MODULE_INSTALL_DIR}\")")
install(CODE
        "EXECUTE_PROCESS ( COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${Truchas_MODULE_DIR}\" \"${Truchas_MODULE_INSTALL_DIR}\")"
)

if (NOT Truchas_CONFIG_INSTALL_DIR)
  set(Truchas_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake/truchas)
endif()

# ---------------------------------------------------------------------------- #
# Source Directories 
# ---------------------------------------------------------------------------- #

# External Projects
add_subdirectory(external)

# Main source directory
add_subdirectory(src)

# Old regression tests
add_subdirectory(regressiontests)

# Old TBrook parser crap
if (ENABLE_TBROOK)
  add_subdirectory(tbrookparser)
endif()

add_subdirectory(test)

# ---------------------------------------------------------------------------- #
# Export Definitions 
# ---------------------------------------------------------------------------- #

# Create the target file
#set(Truchas_TARGETS_FILE truchas-targets.cmake)
#export(truchas truchas_exe FILE ${Truchas_TARGETS_FILE})

# Create truchas-config.cmake (find_package(Truchas))
#set(Truchas_CONFIG_FILE ${Truchas_BINARY_DIR}/truchas-config.cmake)
#configure_file(${Truchas_SOURCE_DIR}/truchas-config.in
#               ${Truchas_CONFIG_FILE}
#	       @ONLY)
#install(FILES ${Truchas_CONFIG_FILE} ${Truchas_TARGETS_FILE} 
#        DESTINATION ${Truchas_CONFIG_INSTALL_DIR})

# ---------------------------------------------------------------------------- #
# Package Definitions 
# ---------------------------------------------------------------------------- #

# This will add a 'package' and 'package_source' targets
include(CPack)

# Package name
set(CPACK_PACKAGE_NAME truchas)

# Version information
set(CPACK_PACKAGE_VERSION_MAJOR "${Truchas_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${Truchas_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${Truchas_PATCH_LEVEL}")
set(CPACK_PACKAGE_VERSION       "${Truchas_VERSION}")

# Basic package information
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Truchas: Cool code package")
set(CPACK_PACKAGE_VENDOR "Los Alamos National Laboratory")
set(CPACK_RESOURCE_FILE_LICENSE ${Truchas_SOURCE_DIR}/LICENSE)

# Gzip tar files are defined as my generator, many other options are available
set(CPACK_SOURCE_PACKAGE_FILE_NAME "truchas-${Truchas_VERSION}"
    CACHE INTERNAL "Distribution base file name")


  
# Files to ignore regular expressions are allowed
set(CPACK_SOURCE_IGNORE_FILES
     ".hg"
     ".hg/*"
     "CVS/"
     ".cvs"
     ".cvs/*"
     ".svn"
     ".svn/*"
     ".*.swp"
     "^${Truchas_SOURCE_DIR}/install"
     "^${Truchas_SOURCE_DIR}/build*"
     "Makefile"
     "CMakeFiles"
     "CMakeCache.txt"
     "cmake_install.cmake")


# ---------------------------------------------------------------------------- #
# Status 
# ---------------------------------------------------------------------------- #

# To print out the flags correctly need the upper case string
if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UC)
endif()

# Remove the ; for a nice print out
foreach(i ${Truchas_PREPROCESSOR})
  set(Truchas_PREPROCESSOR_STRING "${Truchas_PREPROCESSOR_STRING} ${i}")
endforeach()

# Status Report
message(STATUS) 
message(STATUS "-------------------------------------------------------") 
message(STATUS)
message(STATUS "Version: ${Truchas_VERSION}")
message(STATUS)
message(STATUS "Build identifier: ${Truchas_BUILD_ID}")
message(STATUS)
message(STATUS "Build date: ${Truchas_BUILD_DATE}")
message(STATUS)
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS)
if (CMAKE_BUILD_TYPE)
message(STATUS "CMAKE_BUILD_TYPE:\t${CMAKE_BUILD_TYPE}")
message(STATUS)
endif(CMAKE_BUILD_TYPE)
message(STATUS "Compiler Define Flag:\t${Truchas_Compiler_Define}")
message(STATUS)
message(STATUS "C Compiler")
message(STATUS "\t${CMAKE_C_COMPILER}")
message(STATUS "\tID:${CMAKE_C_COMPILER_ID}")
message(STATUS "\tVERSION:${CMAKE_C_COMPILER_VERSION}")
if (CMAKE_BUILD_TYPE_UC)
  message(STATUS "\tCMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC}\t"
                      "=${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
endif()		    
if (CMAKE_C_FLAGS)
  message(STATUS "\tCMAKE_C_FLAGS\t"
                      "=${CMAKE_C_FLAGS}")
endif()	
message(STATUS)
message(STATUS "CXX Compiler")
message(STATUS "\t${CMAKE_CXX_COMPILER}")
message(STATUS "\tID:${CMAKE_CXX_COMPILER_ID}")
message(STATUS "\tVERSION:${CMAKE_CXX_COMPILER_VERSION}")
if (CMAKE_BUILD_TYPE_UC)
  message(STATUS "\tCMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UC}\t"
                      "=${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
endif()		    
if (CMAKE_CXX_FLAGS)
  message(STATUS "\tCMAKE_CXX_FLAGS\t"
                      "=${CMAKE_CXX_FLAGS}")
endif()	
message(STATUS)
message(STATUS "Fortran Compiler")
message(STATUS "\t${CMAKE_Fortran_COMPILER}")
message(STATUS "\tID:${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "\tVERSION:${CMAKE_Fortran_COMPILER_VERSION}")
if (CMAKE_BUILD_TYPE_UC)
  message(STATUS "\tCMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}\t"
                      "=${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
endif()		    
if (CMAKE_Fortran_FLAGS)
  message(STATUS "\tCMAKE_Fortran_FLAGS\t"
                      "=${CMAKE_Fortran_FLAGS}")
endif()		
message(STATUS)
if(ENABLE_MPI)
message(STATUS "MPI")
message(STATUS "\tC include and libraries found: ${MPI_C_FOUND}")
if(MPI_C_FOUND)
message(STATUS "\tC Include path(s):${MPI_C_INCLUDE_PATH}")
message(STATUS "\tC Libraries: ${MPI_C_LIBRARIES}")
endif()
message(STATUS "\tFortran libraries found: ${MPI_Fortran_FOUND}")
if(MPI_C_FOUND)
message(STATUS "\tFortran Libraries: ${MPI_Fortran_LIBRARIES}")
endif()
if(ENABLE_TESTS)
message(STATUS "\tMPIEXEC: ${MPIEXEC}")
message(STATUS "\tMPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
endif()
endif()
message(STATUS)
message(STATUS "Python")
message(STATUS "\tPYTHON_VERSION=${PYTHON_VERSION}")
message(STATUS "\tPYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}")
message(STATUS "\tPYTHON_LIBRARIES=${PYTHON_LIBRARIES}")
message(STATUS)
message(STATUS "NumPy")
message(STATUS "\tNUMPY_VERSION=${NUMPY_VERSION}")
message(STATUS "\tNUMPY_INCLUDE_DIRS=${NUMPY_INCLUDE_DIRS}")
message(STATUS)
message(STATUS "LAPACK")
message(STATUS "\tLAPACK_FOUND=${LAPACK_FOUND}")
message(STATUS "\tLAPACK_LINKER_FLAGS=${LAPACK_LINKER_FLAGS}")
message(STATUS "\tLAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
message(STATUS)
message(STATUS "BLAS")
message(STATUS "\tBLAS_FOUND=${BLAS_FOUND}")
message(STATUS "\tBLAS_LINKER_FLAGS=${BLAS_LINKER_FLAGS}")
message(STATUS "\tBLAS_LIBRARIES=${BLAS_LIBRARIES}")
message(STATUS)
message(STATUS)
message(STATUS "External Projects")
message(STATUS)
message(STATUS "Swig")
if(NOT SWIG_BUILD_TARGET)
  message(STATUS "\tSWIG_FOUND=${SWIG_FOUND}")
  message(STATUS "\tSWIG_VERSION=${SWIG_VERSION}")
  message(STATUS "\tSWIG_EXECUTABLE=${SWIG_EXECUTABLE}")
else()
  message(STATUS "\tAdded to the build")
  message(STATUS "\tBuild target: ${SWIG_BUILD_TARGET}")
  message(STATUS "\tInstall directory: ${TruchasExternal_INSTALL_PREFIX}")
endif()  
message(STATUS)
message(STATUS "HYPRE")
if (NOT HYPRE_BUILD_TARGET)
  message(STATUS "\tHYPRE_FOUND=${HYPRE_FOUND}")
  message(STATUS "\tHYPRE_VERSION=${HYPRE_VERSION}")
  message(STATUS "\tHYPRE_IS_PARALLEL=${HYPRE_IS_PARALLEL}")
  message(STATUS "\tHYPRE_INCLUDE_DIRS=${HYPRE_INCLUDE_DIRS}")
  message(STATUS "\tHYPRE_LIBRARIES=${HYPRE_LIBRARIES}")
else()
  message(STATUS "\tAdded to the build")
  message(STATUS "\tBuild target: ${HYPRE_BUILD_TARGET}")
  message(STATUS "\tInstall directory: ${TruchasExternal_INSTALL_PREFIX}")
endif()
message(STATUS)
message(STATUS "HDF5")
if( NOT HDF5_BUILD_TARGET)
  message(STATUS "\tHDF5_FOUND=${HDF5_FOUND}")
  message(STATUS "\tHDF5_VERSION=${HDF5_VERSION}")
  message(STATUS "\tHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR}")
  message(STATUS "\tHDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
  message(STATUS "\tHDF5_C_LIBRARIES=${HDF5_C_LIBRARIES}")
else()  
  message(STATUS "\tAdded to the build")
  message(STATUS "\tBuild Target: ${HDF5_BUILD_TARGET}")
  message(STATUS "\tInstall directory: ${TruchasExternal_INSTALL_PREFIX}")
endif()  
message(STATUS)
message(STATUS "NETCDF")
if(NOT NETCDF_BUILD_TARGET)
  message(STATUS "\tNETCDF_FOUND=${NETCDF_FOUND}")
  message(STATUS "\tNETCDF_VERSION=${NETCDF_VERSION}")
  message(STATUS "\tNETCDF_INCLUDE_DIR=${NETCDF_INCLUDE_DIR}")
  message(STATUS "\tNETCDF_INCLUDE_DIRS=${NETCDF_INCLUDE_DIRS}")
  message(STATUS "\tNETCDF_Fortran_LIBRARIES=${NETCDF_Fortran_LIBRARIES}")
else()  
  message(STATUS "\tAdded to the build")
  message(STATUS "\tBuild Target: ${NETCDF_BUILD_TARGET}")
  message(STATUS "\tInstall directory: ${TruchasExternal_INSTALL_PREFIX}")
endif()  
message(STATUS)
message(STATUS)
message(STATUS "Preprocessor Command:${Truchas_PREPROCESSOR_STRING}")
message(STATUS "-------------------------------------------------------") 


