if(SEARCH_FOR_NETCDF)
  message(STATUS "Searching for a suitable NetCDF library ...")
  find_package(NetCDF "4.1.3" COMPONENTS Fortran)
  if(NETCDF_FOUND)
    if(NOT NETCDF_HAS_NC4)
      message(STATUS "Found unsuitable NetCDF without required netcdf-4 feature")
      set(NETCDF_FOUND False)
    endif()
  endif()
endif()

if(NETCDF_FOUND)
  list(APPEND projects_found "NetCDF")
else()
  list(APPEND projects_to_build "NetCDF")
  set(NETCDF_VERSION "4.1.3")
  if(BUILD_SHARED_LIBS)
    set(netcdf_shlib_flag "--enable-shared" "--disable-static")
  else()
    set(netcdf_shlib_flag "--enable-static" "--disable-shared")
  endif()

  # Defines to get the name-mangling correct for the Fortran interface
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|GNU")
    set(netcdf_cflags "${cflags} -DpgiFortran")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
    set(netcdf_cflags "${cflags} -DNAGf90Fortran")
  endif()

  # Add include directories to find the correct HDF5 header files
  dir_list_to_includes(includes ${HDF5_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})
  set(netcdf_cflags "${netcdf_cflags} ${includes}")

  # Linker library search directories to find the correct HDF5 libraries
  dir_list_to_link_dirs(netcdf_ldflags ${HDF5_LIBRARY_DIRS})
  
  # This may be required instead -- full set of -L/-l link arguments
  #lib_list_to_link_arg(netcdf_ldflags ${HDF5_C_LIBRARIES} ${MPI_C_LIBRARIES})

  externalproject_add(netcdf
    DEPENDS hdf5
    PREFIX netcdf
    URL ${TARFILE_DIR}/netcdf-${NETCDF_VERSION}.tar.gz
    URL_MD5 ead16cb3b671f767396387dcb3c1a814
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
                      CC=${CMAKE_C_COMPILER}
                      CFLAGS=${netcdf_cflags}
                      FC=${CMAKE_Fortran_COMPILER}
                      FCFLAGS=${fflags}
                      LDFLAGS=${netcdf_ldflags}
                      --enable-netcdf-4
                      --disable-examples
                      --disable-dap
                      --disable-cxx
                      ${netcdf_shlib_flag}
                      --prefix=${CMAKE_INSTALL_PREFIX}
    PATCH_COMMAND patch -p1 < ${TARFILE_DIR}/netcdf-4.1.3.patch
    BUILD_COMMAND "$(MAKE)"
    BUILD_IN_SOURCE 1
    LOG_UPDATE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
  # These FindNetCDF variables are needed to configure exodus.
  set(NETCDF_C_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
  if(BUILD_SHARED_LIBS)
    set(NETCDF_C_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(NETCDF_C_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libnetcdf${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
endif()
