# ############################################################################ #
#
#   DANU  
#     Build Required External Packages
#
# ############################################################################ #

# CMake Modules
include(ExternalProject)

# Danu Modules
include(PrintVariable)
include(BuildLibraryName)
include(BuildAutoconfCflags)

include(ExternalPackageVersions)

#------------------------------------------------------------------------------#
# Set common build compiler flags, build types and directories
#------------------------------------------------------------------------------#

# CMake compiler settings for any package built with CMake
set(Danu_CMAKE_COMPILER_ARGS
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
  -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
  -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
  -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
)

# CMake build settings
set(Danu_CMAKE_BUILD_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
)

# GNU configure static/shared
set(Danu_SHARED_SWITCH --disable-shared)
if(BUILD_SHARED_LIBS)
  set(Danu_SHARED_SWITCH --enable-shared)
endif()

#------------------------------------------------------------------------------#
# Common build definitions
#------------------------------------------------------------------------------#

# All archive files found in the source directory under tpl
set(Danu_TPL_DOWNLOAD_DIR ${Danu_SOURCE_DIR}/tpl)

# All build (binary) directories under the current binary directory
set(Danu_TPL_BINARY_DIR   ${Danu_BINARY_DIR}/tpl-builds)

# Log all activity
set(Danu_TPL_LOG_OPTS
    LOG_DOWNLOAD  1
    LOG_UPDATE    1
    LOG_CONFIGURE 1
    LOG_BUILD     1
    LOG_TEST      1
    LOG_INSTALL   1
)

# Default is to install in the CMAKE_INSTALL_PRFIX area
if (NOT Danu_TPL_INSTALL_PREFIX )
  set(Danu_TPL_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# Create separate targets for each step 
set_property(DIRECTORY PROPERTY
             EP_STEP_TARGETS download patch configure build install test)

# Shared or static libraries flag to use when building library names
set(lib_type STATIC)
if(BUILD_SHARED_LIBS)
  set(lib_type SHARED)
endif()

#------------------------------------------------------------------------------#
# Build: ZLIB 
#------------------------------------------------------------------------------#
set(ZLIB_BUILD_TARGET)
if ( BUILD_ZLIB )

  set(ZLIB_BUILD_TARGET zlib)
  ExternalProject_Add(${ZLIB_BUILD_TARGET}
                      URL     ${Danu_TPL_DOWNLOAD_DIR}/${ZLIB_ARCHIVE_FILE}
                      URL_MD5 ${ZLIB_MD5SUM}
                      UPDATE_COMMAND ""
                      PREFIX  ${Danu_TPL_BINARY_DIR}
                      INSTALL_DIR ${Danu_TPL_INSTALL_PREFIX}
                      CMAKE_CACHE_ARGS
                                ${Danu_CMAKE_COMPILER_ARGS}
                                ${Danu_CMAKE_BUILD_ARGS}
                               -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                      CMAKE_ARGS
                               -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                      ${Danu_TPL_LOG_OPTS})       

  # Define include dirs
  set(ZLIB_INCLUDE_DIRS ${Danu_TPL_INSTALL_PREFIX}/include)

  # Define the zlib libraries
  build_library_name(z ZLIB_LIBRARIES
                     APPEND_PATH ${Danu_TPL_INSTALL_PREFIX}/lib ${lib_type})

else(BUILD_ZLIB)

  if ( BUILD_HDF5 )
    find_package(ZLIB REQUIRED)
  endif()  

endif(BUILD_ZLIB)  

#------------------------------------------------------------------------------#
# Build: HDF5
#------------------------------------------------------------------------------#

set(HDF5_BUILD_TARGET)
if ( BUILD_HDF5 )
  
   find_library(MATH_LIBRARY m)
 
   set(HDF5_BUILD_TARGET hdf5) 

   option(BUILD_HDF5Cmake "Build HDF5 with CMake" OFF)
   
   if ( BUILD_HDF5Cmake ) 
     ExternalProject_Add(${HDF5_BUILD_TARGET}
                         DEPENDS ${ZLIB_BUILD_TARGET}
                         URL ${Danu_TPL_DOWNLOAD_DIR}/${HDF5_ARCHIVE_FILE}
                         URL_MD5 ${HDF5_MD5SUM}
                         PREFIX ${Danu_TPL_BINARY_DIR}
                         INSTALL_DIR ${Danu_TPL_INSTALL_PREFIX}
                         CMAKE_CACHE_ARGS
                                    ${Danu_CMAKE_COMPILER_ARGS}
                                    ${Danu_CMAKE_BUILD_ARGS}
                                   -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                                   -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
                                   -DZIB_INCLUDE_DIR:STRING=${Danu_TPL_INSTALL_PREFIX}/include
                                   -DZLIB_LIBRARY:STRING=${ZLIB_LIBRARIES}
                                   -DHDF5_BUILD_HL_LIB:BOOL=ON
                                   -DHDF5_BUILD_TOOLS:BOOL=ON
                                   -DHDF5_BUILD_FORTRAN:BOOL=OFF
                                   -DHDF5_BUILD_CPP_LIB:BOOL=OFF
                         CMAKE_ARGS
                                   -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                         BUILD_COMMAND $(MAKE)           
                         ${Danu_TPL_LOG_OPTS})
  else()
    build_autoconf_cflags(hdf5_cflags)

    include(BuildWhitespaceString)
    if ( HDF5_LINK_FLAGS )
      build_whitespace_string(hdf5_ldflags -L${Danu_TPL_INSTALL_PREFIX}/lib ${HDF5_LINK_FLAGS})
    else()  
      set(hdf5_ldflags -L${Danu_TPL_INSTALL_PREFIX}/lib)
    endif()

    ExternalProject_Add(${HDF5_BUILD_TARGET}
                         DEPENDS ${ZLIB_BUILD_TARGET}
                         URL ${Danu_TPL_DOWNLOAD_DIR}/${HDF5_ARCHIVE_FILE}
                         URL_MD5 ${HDF5_MD5SUM}
                         PREFIX ${Danu_TPL_BINARY_DIR}
                         INSTALL_DIR ${Danu_TPL_INSTALL_PREFIX}
                         CONFIGURE_COMMAND 
                             <SOURCE_DIR>/configure
                                          --prefix=<INSTALL_DIR>
                                          --enable-option-checking
                                          --disable-fortran
                                          --disable-cxx
                                          --disable-parallel
                                          --enable-largefile
                                          --enable-production
                                          --enable-hl
                                          --with-zlib=${Danu_TPL_INSTALL_PREFIX}
                                          ${Danu_SHARED_SWITCH}
					  CC=${CMAKE_C_COMPILER}
                                          CFLAGS=${hdf5_cflags}
                                          LDFLAGS=${hdf5_ldflags}
					  PARALLEL=
					  RUNPARALLEL=
                          BUILD_COMMAND $(MAKE)                
                          ${Danu_TPL_LOG_OPTS})
  endif()



  # Construct the include variables
  set(HDF5_INCLUDE_DIR ${Danu_TPL_INSTALL_PREFIX}/include)
  set(HDF5_INCLUDE_DIRS ${Danu_TPL_INSTALL_PREFIX}/include)
 
  set(HDF5_LIBRARY_DIRS ${Danu_TPL_INSTALL_PREFIX}/lib)

  # Debug or not, need to worry about this with CMake builds
  set(debug_suffix)
  #string(TOLOWER "${CMAKE_BUILD_TYPE}" build_type_lc) 
  #if ( "${build_type_lc}" STREQUAL "debug"  AND ${BUILD_HDF5Cmake})
  #  set(debug_suffix _debug)
  #endif()

  build_library_name(hdf5_hl${debug_suffix}
                     HDF5_HL_LIBRARY 
                     APPEND_PATH ${Danu_TPL_INSTALL_PREFIX}/lib
                     ${lib_type})

  build_library_name(hdf5${debug_suffix}
                     HDF5_C_LIBRARY
                     APPEND_PATH ${Danu_TPL_INSTALL_PREFIX}/lib
                     ${lib_type})

  set(HDF5_LINK_LIBRARIES ${ZLIB_LIBRARIES} ${MATH_LIBRARY})

  set(HDF5_LIBRARIES 
      ${HDF5_HL_LIBRARY}
      ${HDF5_C_LIBRARY}
      ${HDF5_LINK_LIBRARIES})

endif(BUILD_HDF5)  

#------------------------------------------------------------------------------#
# Build: SWIG
#------------------------------------------------------------------------------#

set(SWIG_BUILD_TARGET)
if ( BUILD_SWIG )

  set(SWIG_BUILD_TARGET swig)

  build_autoconf_cflags(swig_cpp_flags)

  find_package(Python REQUIRED)

  ExternalProject_Add(${SWIG_BUILD_TARGET}
                        URL         ${Danu_TPL_DOWNLOAD_DIR}/${SWIG_ARCHIVE_FILE}
                        URL_MD5     ${SWIG_MD5SUM}
                        PREFIX      ${Danu_TPL_BINARY_DIR}
                        INSTALL_DIR ${Danu_TPL_INSTALL_PREFIX}
                        CONFIGURE_COMMAND
                                    <SOURCE_DIR>/configure
                                       --prefix=<INSTALL_DIR>
                                       --without-pcre
                                       --without-allang
                                       --with-python=${PYTHON_EXECUTABLE}
                                       ${Danu_SHARED_SWITCH}
				       CC=${CMAKE_C_COMPILER}
				       CXX=${CMAKE_CXX_COMPILER}
                                       CPPFLAGS=${swig_cpp_flags}
                        BUILD_COMMAND $(MAKE)           
                        ${Danu_TPL_LOG_OPTS})

  # Define the SWIG related variables
  set(SWIG_DIR          ${Danu_TPL_INSTALL_PREFIX})
  set(SWIG_EXECUTABLE   ${Danu_TPL_INSTALL_PREFIX}/bin/swig)

endif(BUILD_SWIG)


#------------------------------------------------------------------------------#
# Build: Check ( C Unit Test Framework )
#------------------------------------------------------------------------------#

set(Check_BUILD_TARGET)
if ( BUILD_Check )

  set(Check_BUILD_TARGET ucheck)

  build_autoconf_cflags(ucheck_cpp_flags)

  ExternalProject_Add(${Check_BUILD_TARGET}
         URL         ${Danu_TPL_DOWNLOAD_DIR}/${Check_ARCHIVE_FILE}
         URL_MD5     ${Check_MD5SUM}
         PREFIX      ${Danu_TPL_BINARY_DIR}
         INSTALL_DIR ${Danu_TPL_INSTALL_PREFIX}
         CONFIGURE_COMMAND
                    <SOURCE_DIR>/configure
                     --prefix=<INSTALL_DIR>
                     ${Danu_SHARED_SWITCH}
		     CC=${CMAKE_C_COMPILER}
                     CPPFLAGS=${ucheck_cpp_flags}
         BUILD_COMMAND $(MAKE)           
         ${Danu_TPL_LOG_OPTS})

  # Define the include variables
  set(Check_INCLUDE_DIR    ${Danu_TPL_INSTALL_PREFIX}/include)
  set(Check_INCLUDE_DIRS   ${Danu_TPL_INSTALL_PREFIX}/include)
 
  # Define libraries
  build_library_name(check Check_LIBRARIES
                     APPEND_PATH ${Danu_TPL_INSTALL_PREFIX}/lib
                     ${lib_type})

endif(BUILD_Check)  
