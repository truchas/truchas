# #  -*- mode: cmake -*-

#
# ExternalProject_ZLIB
#    Build the ZLIB project
#
#

# --- Setup

# Projects/targets dependent on ZLIB need this target
set(ZLIB_BUILD_TARGET zlib)

# Define the version and archive file
include(ExternalProjectVersions)

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(zlib_url_file     ${TruchasExternal_ARCHIVE_DIR}/${ZLIB_ARCHIVE_FILE})
set(zlib_prefix_dir   ${TruchasExternal_BINARY_DIR}/zlib)
set(zlib_source_dir   ${zlib_prefix_dir}/zlib-${ZLIB_VERSION}-source)
set(zlib_stamp_dir    ${zlib_prefix_dir}/zlib-timestamps)
set(zlib_tmp_dir      ${zlib_prefix_dir}/zlib-tmp)
set(zlib_install_dir  ${TruchasExternal_INSTALL_PREFIX})

# --- Add the -fPIC or -PIC (Position in code flag)
include(FindPICFlag)
find_pic_flag(zlib_pic_flag)
set(zlib_c_flags "${CMAKE_C_FLAGS} ${zlib_pic_flag}")


# --- Add the external project

ExternalProject_Add(${ZLIB_BUILD_TARGET}
                    # -- Project directories
                    PREFIX      ${zlib_prefix_dir}   
                    TMP_DIR     ${zlib_tmp_dir}     
                    STAMP_DIR   ${zlib_stamp_dir}
                    SOURCE_DIR  ${zlib_source_dir}
		    INSTALL_DIR ${zlib_install_dir}
		    # -- Archive file definitions
                    URL          ${zlib_url_file}
                    URL_MD5      ${ZLIB_MD5_SUM}   
                    # -- Configure (CMake)
		    CMAKE_CACHE_ARGS
		          ${TruchasExternal_CMAKE_COMPILER_ARGS}
			  ${TruchasExternal_CMAKE_BUILD_ARGS}
			  -DCMAKE_C_FLAGS:STRING=${zlib_c_flags}
	            CMAKE_ARGS
		          -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                    # -- Output control
                    ${zlib_logging})

# --- Set the variables for other targets that need ZLIB

# Version
global_set(ZLIB_VERSION ${ZLIB_VERSION})

# Include directory
global_set(ZLIB_INCLUDE_DIR ${zlib_install_dir}/include)
global_set(ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS})

# Library
build_library_name(z library APPEND_PATH ${zlib_install_dir}/lib)
global_set(ZLIB_LIBRARY ${library})
global_set(ZLIB_LIBRARIES ${library})
