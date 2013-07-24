# #  -*- mode: cmake -*-

#
# ExternalProject_SWIG
#    Build the SWIG project
#
#

# --- Setup

# Projects/targets dependent on SWIG need this target
set(SWIG_BUILD_TARGET swig)
global_set(SWIG_BUILD_TARGET ${SWIG_BUILD_TARGET})

# Define the version and archive file
include(ExternalProjectVersions)

# Useful utility to build *FLAGS strings
include(BuildWhitespaceString)

# Make CMake files to run scripts and exit correctly 
include(MakeCMakeCommandFile)

# Construct library names
include(BuildLibraryName)

# ExternalProject directories, file and log settings
set(swig_url_file     ${TruchasExternal_ARCHIVE_DIR}/${SWIG_ARCHIVE_FILE})
set(swig_prefix_dir   ${TruchasExternal_BINARY_DIR}/swig)
set(swig_source_dir   ${swig_prefix_dir}/swig-${SWIG_VERSION}-source)
set(swig_stamp_dir    ${swig_prefix_dir}/swig-timestamps)
set(swig_tmp_dir      ${swig_prefix_dir}/swig-tmp)
set(swig_install_dir  ${TruchasExternal_INSTALL_PREFIX})

# --- Configure flags
if(ENABLE_MPI)
  set(swig_parallel_opt --enable-parallel)
else()
  set(swig_parallel_opt --disable-parallel)
endif()

# --- Compile flags
string(TOUPPER "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_FLAGS)

build_whitespace_string(swig_cflags 
                        ${CMAKE_C_FLAGS}
			${${CMAKE_BUILD_TYPE_FLAGS}})

# --- Add the external project

ExternalProject_Add(${SWIG_BUILD_TARGET}
                    # -- Project directories
                    PREFIX      ${swig_prefix_dir}   
                    TMP_DIR     ${swig_tmp_dir}     
                    STAMP_DIR   ${swig_stamp_dir}
                    SOURCE_DIR  ${swig_source_dir}
		    #INSTALL_DIR ${swig_install_dir}
		    # -- Archive file definitions
                    URL          ${swig_url_file}
                    URL_MD5      ${SWIG_MD5_SUM}   
                    # -- Configure
		    CONFIGURE_COMMAND
                        <SOURCE_DIR>/configure
                                          --prefix=${swig_install_dir}
					  --without-pcre
					  --with-python=${PYTHON_EXECUTABLE}
                                          ${TruchasExternal_SHARED_SWITCH}
					  CC=${CMAKE_C_COMPILER}
					  CXX=${CMAKE_CXX_COMPILER}
                                          CPPFLAGS=${swig_cflags}
                    # -- Build					  
                    BUILD_COMMAND $(MAKE)             
                    # -- Output control
                    ${TruchasExternal_LOG_OPTS})

# --- Set the variables for other targets that need SWIG

# Version
global_set(SWIG_VERSION ${SWIG_VERSION})

# Installation directory
global_set(SWIG_DIR ${swig_install_dir})

# SWIG executable
global_set(SWIG_EXECUTABLE ${swig_install_dir}/bin/swig)
