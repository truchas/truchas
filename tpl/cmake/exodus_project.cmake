if(SEARCH_FOR_EXODUS)
  message(STATUS "Searching for a suitable Exodus library ...")
  find_package(Exodus "514")
endif()

if(EXODUS_FOUND)
  list(APPEND projects_found "Exodus")
else()
  list(APPEND projects_to_build "Exodus")
  set(EXODUS_VERSION "514")
  # Go ahead and specify the library paths below so that the find_library calls
  # in the exodus CMakeLists.txt will be short-circuited.  I do not think the
  # libraries actually matter (for shared at least) because we are not building
  # any executables.  Note that an ldd on the built library may not show the
  # correct netcdf and hdf5 libraries, because the ones passed are not baked
  # into it, and whatever is found in LD_LIBRARY_PATH, if any, will be shown.
  # But when the truchas executable is linked, the proper libraries will be
  # used.  The libraries we want are the initial ones in the lists.
  list(GET HDF5_C_LIBRARIES 0 hdf5_library)
  list(GET HDF5_HL_LIBRARIES 0 hdf5hl_library)
  ExternalProject_Add(exodus
    DEPENDS hdf5 netcdf
    PREFIX exodus
    URL ${TARFILE_DIR}/exodus-5.14.tar.gz
    URL_MD5 064ec8ea279b4b7f63bd113645dc3501
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
               -D CMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
               #-D CMAKE_PREFIX_PATH:PATH=${CMAKE_INSTALL_PREFIX}
               -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
               -D NETCDF_INCLUDE_DIR:PATH=${NETCDF_C_INCLUDE_DIR}
               -D NETCDF_LIBRARY:PATH=${NETCDF_C_LIBRARY}
               -D HDF5_LIBRARY=${hdf5_library}
               -D HDF5HL_LIBRARY=${hdf5hl_library}
               #-D CMAKE_EXE_LINKER_FLAGS=${hdf5_hl_ldflags}
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

