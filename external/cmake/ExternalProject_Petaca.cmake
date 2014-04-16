# #  -*- mode: cmake -*-

#
# ExternalProject_Petaca
#    Build the Petaca project
#
#

# --- Setup

# Projects/targets dependent on Petaca need this target
set(Petaca_BUILD_TARGET petaca)

# Define the version and archive file
set(EP_Petaca_VERSION 0768bc4)
set(EP_Petaca_ARCHIVE_FILE   petaca-${EP_Petaca_VERSION}.tar.gz)
set(EP_Petaca_MD5_SUM        3b3088f7f54fdffd9c2eef45436a88b3)  

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(petaca_url_file     ${TruchasExternal_ARCHIVE_DIR}/${EP_Petaca_ARCHIVE_FILE})
set(petaca_prefix_dir   ${TruchasExternal_BINARY_DIR}/petaca)
set(petaca_source_dir   ${petaca_prefix_dir}/petaca-${EP_Petaca_VERSION}-source)
set(petaca_stamp_dir    ${petaca_prefix_dir}/petaca-timestamps)
set(petaca_tmp_dir      ${petaca_prefix_dir}/petaca-tmp)
set(petaca_install_dir  ${TruchasExternal_INSTALL_PREFIX})

# YAJL (Required)
if ( NOT TARGET ${YAJL_BUILD_TARGET} )
  include(Verify_YAJL)
  if (NOT YAJL_VERIFIED)
    include(ExternalProject_YAJL)
  endif(NOT YAJL_VERIFIED)  
endif(NOT TARGET ${YAJL_BUILD_TARGET})

# --- Add the -fPIC or -PIC (Position in code flag)
include(FindPICFlag)
find_pic_flag(petaca_pic_flag)
set(petaca_c_flags "${CMAKE_C_FLAGS} ${petaca_pic_flag}")

# --- Add the external project
ExternalProject_Add(${Petaca_BUILD_TARGET}
                    DEPENDS ${YAJL_BUILD_TARGET}
                    # -- Project directories
                    PREFIX      ${petaca_prefix_dir}   
                    TMP_DIR     ${petaca_tmp_dir}     
                    STAMP_DIR   ${petaca_stamp_dir}
                    SOURCE_DIR  ${petaca_source_dir}
		    #INSTALL_DIR ${petaca_install_dir}
		    # -- Archive file definitions
                    URL          ${petaca_url_file}
                    URL_MD5      ${EP_Petaca_MD5_SUM}   
                    # -- Configure (CMake)
		    CMAKE_CACHE_ARGS
		          ${TruchasExternal_CMAKE_COMPILER_ARGS}
			  ${TruchasExternal_CMAKE_BUILD_ARGS}
			  -DCMAKE_C_FLAGS:STRING=${petaca_c_flags}
	            CMAKE_ARGS
		          -DCMAKE_INSTALL_PREFIX:PATH=${petaca_install_dir}
                          #-DYAJL_INCLUDE_DIR:PATH=${YAJL_INCLUDE_DIR}
                          #-DYAJL_LIBRARY_DIR:PATH=${YAJL_LIBRARY_DIR}
                    # -- Output control
                    ${TruchasExternal_LOG_OPTS})

# --- Set the variables for other targets that need Petaca

# Version
set(Petaca_VERSION ${EP_Petaca_VERSION})

# Include directory
set(Petaca_INCLUDE_DIR ${petaca_install_dir}/include)
set(Petaca_INCLUDE_DIRS ${Petaca_INCLUDE_DIR})

# Library
build_library_name(petaca Petaca_LIBRARY APPEND_PATH ${petaca_install_dir}/lib)
set(Petaca_LIBRARIES ${Petaca_LIBRARY} yajl)

