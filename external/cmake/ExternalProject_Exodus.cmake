# #  -*- mode: cmake -*-

#
# ExternalProject_EXODUS
#    Build the ExodusII library
#
#

# --- Setup

# Project build target
set(EXODUS_BUILD_TARGET exodus)

# Define the version and archive file
set(EP_EXODUS_VERSION 514)
set(EP_EXODUS_ARCHIVE_FILE exodus-5.14.tar.gz)
set(EP_EXODUS_MD5_SUM 064ec8ea279b4b7f63bd113645dc3501) 

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(exodus_url_file     ${TruchasExternal_ARCHIVE_DIR}/${EP_EXODUS_ARCHIVE_FILE})
set(exodus_prefix_dir   ${TruchasExternal_BINARY_DIR}/exodus)
set(exodus_source_dir   ${exodus_prefix_dir}/exodus-source)
set(exodus_binary_dir   ${exodus_prefix_dir}/exodus-build)
set(exodus_stamp_dir    ${exodus_prefix_dir}/exodus-timestamps)
set(exodus_tmp_dir      ${exodus_prefix_dir}/exodus-tmp)
set(exodus_install_dir  ${TruchasExternal_INSTALL_PREFIX})
set(exodus_download_dir ${TruchasExternal_ARCHIVE_DIR})

# --- Add the -fPIC or -PIC (Position in code flag)
include(FindPICFlag)
find_pic_flag(exodus_pic_flag)
set(exodus_c_flags "${CMAKE_C_FLAGS} ${exodus_pic_flag}")

# --- Add the external project

ExternalProject_Add(${EXODUS_BUILD_TARGET}
                    DEPENDS ${NETCDF_BUILD_TARGET} ${HDF5_BUILD_TARGET}
                    # -- Project directories
                    PREFIX      ${exodus_prefix_dir} 
                    TMP_DIR     ${exodus_tmp_dir}     
                    STAMP_DIR   ${exodus_stamp_dir}
                    SOURCE_DIR  ${exodus_source_dir}
                    BINARY_DIR  ${exodus_binary_dir}
                    #INSTALL_DIR ${exodus_install_dir}
                    # -- Archive file definitions
                    DOWNLOAD_DIR ${exodus_download_dir}
                    URL          ${exodus_url_file}
                    URL_MD5      ${EP_EXODUS_MD5_SUM}   
                    # -- Configure (CMake)
                    CMAKE_COMMAND ${CMAKE_COMMAND} ${exodus_source_dir}/exodus
                    CMAKE_CACHE_ARGS
                        ${TruchasExternal_CMAKE_COMPILER_ARGS}
                        ${TruchasExternal_CMAKE_BUILD_ARGS}
                       -DCMAKE_C_FLAGS:STRING=${exodus_c_flags}
                    CMAKE_ARGS
                       -DCMAKE_INSTALL_PREFIX:PATH=${exodus_install_dir}
                    # -- Output control
                    ${Truchas_External_LOG_OPTS})

# --- Set the variables for other targets that need exodus

# Version
set(EXODUS_VERSION ${EP_EXODUS_VERSION})

# Include directory
set(EXODUS_INCLUDE_DIR ${exodus_install_dir}/include)
set(EXODUS_INCLUDE_DIRS ${EXODUS_INCLUDE_DIR})

# Library
build_library_name(exoIIv2c EXODUS_LIBRARY APPEND_PATH ${exodus_install_dir}/lib)
set(EXODUS_LIBRARIES ${EXODUS_LIBRARY} ${NETCDF_C_LIBRARIES})
