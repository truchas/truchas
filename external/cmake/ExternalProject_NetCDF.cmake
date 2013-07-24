# #  -*- mode: cmake -*-

#
# ExternalProject_NetCDF
#    Build the NetCDF project
#
#

# --- Setup

# Projects/targets dependent on NetCDF need this target
set(NetCDF_BUILD_TARGET netcdf)
global_set(NetCDF_BUILD_TARGET ${NetCDF_BUILD_TARGET})

# Define the version and archive file
include(ExternalProjectVersions)

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(netcdf_url_file     ${TruchasExternal_ARCHIVE_DIR}/${NetCDF_ARCHIVE_FILE})
set(netcdf_prefix_dir   ${TruchasExternal_BINARY_DIR}/netcdf)
set(netcdf_source_dir   ${netcdf_prefix_dir}/netcdf-${NetCDF_VERSION}-source)
set(netcdf_stamp_dir    ${netcdf_prefix_dir}/netcdf-timestamps)
set(netcdf_tmp_dir      ${netcdf_prefix_dir}/netcdf-tmp)
set(netcdf_install_dir  ${TruchasExternal_INSTALL_PREFIX})

# Need the Compiler names without the full path
get_filename_component(CMAKE_C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component(CMAKE_Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# --- Patch defines

# Need Perl
find_package(Perl)
if (NOT PERL_FOUND)
  message(FATAL_ERROR "Can not locate Perl. Can not patch NetCDF.")
endif()

# Configure the bash patch script
set(netcdf_sh_patch ${netcdf_prefix_dir}/netcdf-patch-step.sh)
configure_file(${TruchasExternal_TEMPLATES_DIR}/netcdf-patch-step.sh.in
               ${netcdf_sh_patch}
               @ONLY)

# Configure the  CMake command file
set(netcdf_cmake_patch ${netcdf_prefix_dir}/netcdf-patch-step.cmake)
make_cmake_command_file(${netcdf_cmake_patch}
                        SCRIPT ${netcdf_sh_patch}
                        ACTION patch
			EP_NAME NetCDF
			WORK_DIRECTORY ${netcdf_prefix_dir})

set(NetCDF_PATCH_COMMAND ${CMAKE_COMMAND} -P ${netcdf_cmake_patch})  

# --- Compile flags
string(TOUPPER CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} CMAKE_C_BUILD_FLAGS)
string(TOUPPER CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} CMAKE_CXX_BUILD_FLAGS)
string(TOUPPER CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} CMAKE_Fortran_BUILD_FLAGS)

# C flags
if ( "${CMAKE_Fortran_COMPILER_NAME}" MATCHES "Nag|nag|NAG" )
  set(fortran_define -DNAGf90Fortran)
elseif ( "${CMAKE_Fortran_COMPILER_NAME}" MATCHES "g95" )
  set(fortran_define -DNAGf90Fortran)
endif()
build_whitespace_string(netcdf_cflags 
                        ${CMAKE_C_FLAGS}
			${${CMAKE_C_BUILD_FLAGS}}
			${fortran_define})

# CPP flags
build_whitespace_string(netcdf_cppflags 
                        -O -DNDEBUG)         

# Fortran flags
#build_whitespace_string(netcdf_fcflags 
#                        ${CMAKE_Fortran_FLAGS}
#			${${CMAKE_Fortran_BUILD_FLAGS}}
#			)
build_whitespace_string(netcdf_fcflags 
                        -O -DNDEBUG)


# --- Link flags
build_whitespace_string(netcdf_ldflags -L${zlib_install_dir}/lib)


# --- Add the external project

ExternalProject_Add(${NetCDF_BUILD_TARGET}
                    DEPENDS ${ZLIB_BUILD_TARGET}
                    # -- Project directories
                    PREFIX      ${netcdf_prefix_dir}   
                    TMP_DIR     ${netcdf_tmp_dir}     
                    STAMP_DIR   ${netcdf_stamp_dir}
                    SOURCE_DIR  ${netcdf_source_dir}
		    #INSTALL_DIR ${netcdf_install_dir}
		    # -- Archive file definitions
                    URL          ${netcdf_url_file}
                    URL_MD5      ${NetCDF_MD5_SUM}   
		    # -- Patch
		    PATCH_COMMAND ${NetCDF_PATCH_COMMAND}
                    # -- Configure
		    CONFIGURE_COMMAND <SOURCE_DIR>/configure
		                          --prefix=${netcdf_install_dir}
		        		  --disable-dap
		    			  --disable-netcdf-4
		    			  --disable-cxx
		    			  --with-pic
                                          ${TruchasExternal_SHARED_SWITCH}
		    			  CC=${CMAKE_C_COMPILER_NAME}
		    			  CXX=""
		    			  FC=${CMAKE_Fortran_COMPILER_NAME}
		                         CFLAGS=${netcdf_cflags}
		                          CPPFLAGS=${netcdf_cppflags}
		    			  FCFLAGS=${netcdf_fcflags}
		                         LDFLAGS=${netcdf_ldflags}
                    # -- Build					  
                    BUILD_COMMAND $(MAKE)             
		    BUILD_IN_SOURCE True
                    # -- Output control
                    ${TruchasExternal_LOG_OPTS})

# --- Set the variables for other targets that need NetCDF

# Version
global_set(NetCDF_VERSION ${NetCDF_VERSION})

# Include directory
global_set(NetCDF_INCLUDE_DIR ${netcdf_install_dir}/include)
global_set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})

# Library
global_set(NetCDF_LIBRARY_DIRS ${netcdf_install_dir}/lib)

build_library_name(netcdf
                   NetCDF_C_LIBRARY 
                   APPEND_PATH ${netcdf_install_dir}/lib)
global_set(NetCDF_C_LIBRARY ${NetCDF_C_LIBRARY})

build_library_name(netcdff
                   NetCDF_Fortran_LIBRARY 
                   APPEND_PATH ${netcdf_install_dir}/lib)
global_set(NetCDF_Fortran_LIBRARY ${NetCDF_Fortran_LIBRARY})

global_set(NetCDF_C_LIBRARIES 
           ${NetCDF_C_LIBRARY} ${ZLIB_LIBRARIES})

global_set(NetCDF_Fortran_LIBRARIES 
           ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARIES})


