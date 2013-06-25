# #  -*- mode: cmake -*-

#
# ExternalProject_HYPRE
#    Build the HYPRE project
#
#

# --- Setup

# Projects/targets dependent on HYPRE need this target
set(HYPRE_BUILD_TARGET hypre)

# Define the version and archive file
include(ExternalProjectVersions)

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(hypre_url_file     ${TruchasExternal_ARCHIVE_DIR}/${HYPRE_ARCHIVE_FILE})
if (ENABLE_MPI)
  set(hypre_prefix_dir   ${TruchasExternal_BINARY_DIR}/hypre-mpi)
else()
  set(hypre_prefix_dir   ${TruchasExternal_BINARY_DIR}/hypre-serial)
endif()  
set(hypre_source_dir   ${hypre_prefix_dir}/hypre-${HYPRE_VERSION}-source)
set(hypre_stamp_dir    ${hypre_prefix_dir}/hypre-timestamps)
set(hypre_tmp_dir      ${hypre_prefix_dir}/hypre-tmp)
set(hypre_install_dir  ${TruchasExternal_INSTALL_PREFIX})

# --- Configure parameters

# MPI flags and opts
if ( ENABLE_MPI )
  build_whitespace_string(mpi_include ${MPI_C_INCLUDE_PATH})
  build_whitespace_string(mpi_libs ${MPI_C_LIBRARIES})
  build_whitespace_string(hypre_mpi_opt
                          --with-MPI 
                          --with-MPI-include='${mpi_include}'
			  --with-MPI-libs='${mpi_libs}')
  set(hypre_mpi_cflags -I${MPI_C_INCLUDE_PATH})
else()
  set(hypre_mpi_opt --without-MPI)
  set(hypre_mpi_cflags)
endif()

#BROKEN# BLAS options
#BROKENset(hypre_blas_opt)
#BROKENfind_package(BLAS)
#BROKENif (BLAS_FOUND)
#BROKEN  build_whitespace_string(blas_libs ${BLAS_LIBRARIES}) 
#BROKEN  build_whitespace_string(hypre_blas_opt
#BROKEN                         --with-blas-libs='${blas_libs}')
#BROKENelse()  
#BROKENendif()
#BROKEN
#BROKEN# LAPACK options
#BROKENset(hypre_lapack_opt)
#BROKENfind_package(LAPACK)
#BROKENif (LAPACK_FOUND)
#BROKEN  build_whitespace_string(lapack_libs ${BLAS_LIBRARIES}) 
#BROKEN  build_whitespace_string(hypre_lapack_opt
#BROKEN                         --with-lapack-libs='${blas_libs}')
#BROKENendif()
#BROKEN
#BROKEN# OpenMP
#BROKENset(hypre_openmp_opt)
#BROKENset(hypre_openmp_flags)
#BROKENfind_package(OpenMP)
#BROKENif (OPENMP_FOUND)
#BROKEN  set(hypre_openmp_opt --with-openmp)
#BROKEN  set(hypre_openmp_flags ${OpenMP_C_FLAGS})
#BROKENendif()

# Compile flags 
build_whitespace_string(hypre_cflags
                        ${hypre_mpi_cflags}
			${hypre_openmp_flags})

# --- Create the configure scripts and command		      

# Build the configure script
set(hypre_sh_configure ${hypre_prefix_dir}/hypre-configure-step.sh)
configure_file(${TruchasExternal_SOURCE_DIR}/templates/hypre-configure-step.sh.in
               ${hypre_sh_configure}
	       @ONLY)

# Configure the CMake command file
set(hypre_cmake_configure ${hypre_prefix_dir}/hypre-configure-step.cmake)
make_cmake_command_file(${hypre_cmake_configure}
                        SCRIPT ${hypre_sh_configure}
                        ACTION configure
                        EP_NAME HYPRE
		        WORK_DIRECTORY ${hypre_source_dir}/src)

# Configure command	     
set(HYPRE_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${hypre_cmake_configure})	
		      
# --- Define the build scripts and command

# Build script
set(hypre_sh_build ${hypre_prefix_dir}/hypre-build-step.sh)
configure_file(${TruchasExternal_SOURCE_DIR}/templates/hypre-build-step.sh.in
               ${hypre_sh_build}
	       @ONLY)

# CMake build command file
set(hypre_cmake_build ${hypre_prefix_dir}/hypre-build-step.cmake)
make_cmake_command_file(${hypre_cmake_build}
                        SCRIPT ${hypre_sh_build}
                        ACTION build
			EP_NAME HYPRE
			WORK_DIRECTORY ${hypre_source_dir}/src)

# Build command     
set(HYPRE_BUILD_COMMAND ${CMAKE_COMMAND} -P ${hypre_cmake_build})	

# --- Define the install scripts and command

# Build the install script
set(hypre_sh_install ${hypre_prefix_dir}/hypre-install-step.sh)
configure_file(${TruchasExternal_SOURCE_DIR}/templates/hypre-install-step.sh.in
               ${hypre_sh_install}
	       @ONLY)

# CMake install command file
set(hypre_cmake_install ${hypre_prefix_dir}/hypre-install-step.cmake)
make_cmake_command_file(${hypre_cmake_install}
                        SCRIPT ${hypre_sh_install}
                        ACTION install
			EP_NAME HYPRE
			WORK_DIRECTORY ${hypre_source_dir}/src)

# Install command
set(HYPRE_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${hypre_cmake_install})	

# --- Add the external project
ExternalProject_Add(${HYPRE_BUILD_TARGET}
                    # -- Project directories
                    PREFIX     ${hypre_prefix_dir}   
                    TMP_DIR    ${hypre_tmp_dir}     
                    STAMP_DIR  ${hypre_stamp_dir}
		    # -- Archive file definitions
                    URL          ${hypre_url_file}
                    URL_MD5      ${HYPRE_MD5_SUM}   
                    # -- Configure
                    SOURCE_DIR        ${hypre_source_dir}
                    CONFIGURE_COMMAND ${HYPRE_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${hypre_build_dir}        
                    BUILD_COMMAND     ${HYPRE_BUILD_COMMAND}   
                    BUILD_IN_SOURCE   TRUE                    
                    # -- Install
                    INSTALL_DIR      ${hypre_install_dir}
		    INSTALL_COMMAND  ${HYPRE_INSTALL_COMMAND}
                    # -- Output control
		    ${TruchasExternal_LOG_OPTS})

# --- Set the variables for other targets that need HYPRE

# Version
global_set(HYPRE_VERSION ${HYPRE_VERSION})

# Include directory
global_set(HYPRE_INCLUDE_DIR ${hypre_install_dir}/include)
set(inc_dirs ${HYPRE_INCLUDE_DIR})
if (ENABLE_MPI)
  list(APPEND inc_dir ${MPI_C_INCLUDE_PATH})
endif()
global_set(HYPRE_INCLUDE_DIRS ${inc_dirs})

# Library
build_library_name(HYPRE HYPRE_LIBRARY APPEND_PATH ${hypre_install_dir}/lib)
global_set(HYPRE_LIBRARY ${HYPRE_LIBRARY})
set(libs ${HYPRE_LIBRARY})
if(ENABLE_MPI)
  list(APPEND libs ${MPI_C_LIBRARIES})
endif()
global_set(HYPRE_LIBRARIES ${libs})

# Flags
if (ENABLE_MPI)
  global_set(HYPRE_IS_PARALLEL True)
else()
  global_set(HYPRE_IS_PARALLEL False)
endif()

