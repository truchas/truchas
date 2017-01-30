if(SEARCH_FOR_HYPRE)
  message(STATUS "Searching for a suitable HYPRE library ...")
  find_package(HYPRE 2.6.0 EXACT)
  if(HYPRE_FOUND)
    if(NOT HYPRE_IS_PARALLEL)
      set(HYPRE_FOUND False)
      message(STATUS "Require parallel HYPRE library but found unsuitable serial library")
    endif()
  endif()
endif()

if(HYPRE_FOUND)
  list(APPEND projects_found "HYPRE")
else()
  list(APPEND projects_to_build "HYPRE")
  set(HYPRE_VERSION "2.6.0b")
  if(BUILD_SHARED_LIBS)
    set(hypre_shlib_flag "--enable-shared")
  else()
    set(hypre_shlib_flag "--disable-shared")
  endif()
  externalproject_add(hypre
    PREFIX hypre
    URL ${TARFILE_DIR}/hypre-${HYPRE_VERSION}.tar.gz
    URL_MD5 84381005bdddff69b62b43ca025070fd
    UPDATE_COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${PROJECT_BINARY_DIR}/hypre/src/hypre/src
        ${PROJECT_BINARY_DIR}/hypre/src/hypre
    CONFIGURE_COMMAND ./configure
                      CC=${MPI_C_COMPILER}
                      CFLAGS=${cflags}
                      --with-MPI
                      --without-fei
                      --disable-fortran
                      --prefix=${CMAKE_INSTALL_PREFIX}
                      ${hypre_shlib_flag}
    BUILD_COMMAND "$(MAKE)" all
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()
