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
set(EP_HYPRE_VERSION_MAJOR  2)
set(EP_HYPRE_VERSION_MINOR  6)
set(EP_HYPRE_VERSION_PATCH  0b)
set(EP_HYPRE_VERSION  ${EP_HYPRE_VERSION_MAJOR}.${EP_HYPRE_VERSION_MINOR}.${EP_HYPRE_VERSION_PATCH})
set(EP_HYPRE_ARCHIVE_FILE   hypre-${EP_HYPRE_VERSION}.tar.gz)
set(EP_HYPRE_MD5_SUM        84381005bdddff69b62b43ca025070fd) 

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(hypre_url_file     ${TruchasExternal_ARCHIVE_DIR}/${EP_HYPRE_ARCHIVE_FILE})
if (ENABLE_MPI)
  set(hypre_prefix_dir   ${TruchasExternal_BINARY_DIR}/hypre-mpi)
else()
  set(hypre_prefix_dir   ${TruchasExternal_BINARY_DIR}/hypre-serial)
endif()  
set(hypre_source_dir   ${hypre_prefix_dir}/hypre-${EP_HYPRE_VERSION}-source)
set(hypre_stamp_dir    ${hypre_prefix_dir}/hypre-timestamps)
set(hypre_tmp_dir      ${hypre_prefix_dir}/hypre-tmp)
set(hypre_install_dir  ${TruchasExternal_INSTALL_PREFIX})
set(hypre_download_dir ${TruchasExternal_ARCHIVE_DIR})

# --- Compile/Build Environment Variables
set(cflags_list)
set(ldflags_list)

# MPI flags and opts
if ( ENABLE_MPI )

  # Configure option
  set(hypre_mpi_opt --with-MPI)

  # MPI flags for CFLAGS
  foreach(dir ${MPI_C_INCLUDE_PATH})
    list(APPEND cflags_list -I${dir})
  endforeach()
  list(APPEND cflags_list ${MPI_C_COMPILE_FLAGS})

  # MPI flags for LDFLAGS
  include(SplitLibraryName)
  foreach(lib ${MPI_C_LIBRARIES})
    split_library_name(${lib} PATH lib_path LIBNAME lib_name EXT lib_ext)
    list(APPEND ldflags_list -L${lib_path} -l${lib_name})
  endforeach()
  list(APPEND ldflags_list ${MPI_C_LINK_FLAGS})

else()
  set(hypre_mpi_opt --without-MPI)
endif()

# BLAS options
set(hypre_blas_opt)
find_package(BLAS)
if (BLAS_FOUND)
  list(APPEND ldflags_list ${BLAS_LIBRARIES})
  list(APPEND ldflags_list ${BLAS_LINKER_FLAGS})
  set(hypre_blas_opt --with-blas)
endif()  

# LAPACK options
set(hypre_lapack_opt)
find_package(LAPACK)
if (LAPACK_FOUND)
  list(APPEND ldflags_list ${LAPACK_LIBRARIES})
  list(APPEND ldflags_list ${LAPACK_LINKER_FLAGS})
  set(hypre_lapack_opt --with-lapack)
endif()

# OpenMP
#BROKENset(hypre_openmp_opt)
#BROKENfind_package(OpenMP)
#BROKENif (OPENMP_FOUND)
#BROKEN  set(hypre_openmp_opt --with-openmp)
#BROKEN  set(cflags_list ${OpenMP_C_FLAGS})
#BROKENendif()

# Whitespace strings for the shell scripts
build_whitespace_string(hypre_cflags ${cflags_list})
build_whitespace_string(hypre_ldflags ${ldflags_list})

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
                    DOWNLOAD_DIR ${hypre_download_dir}
                    URL          ${hypre_url_file}
                    URL_MD5      ${EP_HYPRE_MD5_SUM}   
                    # -- Configure
                    SOURCE_DIR        ${hypre_source_dir}
                    CONFIGURE_COMMAND ${HYPRE_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${hypre_build_dir}        
                    BUILD_COMMAND     ${HYPRE_BUILD_COMMAND}   
                    BUILD_IN_SOURCE   TRUE                    
                    # -- Install
		                #INSTALL_DIR      ${hypre_install_dir}
		                INSTALL_COMMAND  ${HYPRE_INSTALL_COMMAND}
                    # -- Output control
		    ${TruchasExternal_LOG_OPTS})

# --- Set the variables for other targets that need HYPRE

# Version
set(HYPRE_VERSION ${EP_HYPRE_VERSION})

# Include directory
set(HYPRE_INCLUDE_DIR ${hypre_install_dir}/include)
set(inc_dirs ${HYPRE_INCLUDE_DIR})
if (ENABLE_MPI)
  list(APPEND inc_dir ${MPI_C_INCLUDE_PATH})
endif()
set(HYPRE_INCLUDE_DIRS ${inc_dirs})

# Library
build_library_name(HYPRE HYPRE_LIBRARY APPEND_PATH ${hypre_install_dir}/lib)
set(HYPRE_LIBRARY ${HYPRE_LIBRARY})
set(libs ${HYPRE_LIBRARY})
if(ENABLE_MPI)
  list(APPEND libs ${MPI_C_LIBRARIES})
endif()
if(LAPACK_FOUND)
  list(APPEND libs ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
  list(APPEND libs ${BLAS_LIBRARIES})
endif()
set(HYPRE_LIBRARIES ${libs})

# Flags
if (ENABLE_MPI)
  set(HYPRE_IS_PARALLEL True)
else()
  set(HYPRE_IS_PARALLEL False)
endif()

