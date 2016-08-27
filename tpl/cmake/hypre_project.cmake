if(SEARCH_FOR_HYPRE)
  message(STATUS "Searching for a suitable HYPRE library ...")
  find_package(HYPRE 2.6.0 EXACT)
  if(HYPRE_FOUND)
    if(ENABLE_MPI)
      if(NOT HYPRE_IS_PARALLEL)
        set(HYPRE_FOUND False)
        message(STATUS "Require parallel HYPRE library but found unsuitable serial library")
      endif()
    else()
      if(HYPRE_IS_PARALLEL)
        set(HYPRE_FOUND False)
        message(STATUS "Require serial HYPRE library but found unsuitable parallel library")
      endif()
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
  if(ENABLE_MPI)
    set(hypre_mpi_flag "--with-MPI")
    set(hypre_c_compiler ${MPI_C_COMPILER})
  else()
    set(hypre_mpi_flag "--without-MPI")
    set(hypre_c_compiler ${CMAKE_C_COMPILER})
  endif()
  externalproject_add(hypre
    PREFIX hypre
    URL ${TARFILE_DIR}/hypre-${HYPRE_VERSION}.tar.gz
    URL_MD5 84381005bdddff69b62b43ca025070fd
    UPDATE_COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${PROJECT_BINARY_DIR}/hypre/src/hypre/src
        ${PROJECT_BINARY_DIR}/hypre/src/hypre
    CONFIGURE_COMMAND ./configure
                      CC=${hypre_c_compiler}
                      CFLAGS=${cflags}
                      --without-fei
                      --disable-fortran
                      --prefix=${CMAKE_INSTALL_PREFIX}
                      ${hypre_shlib_flag}
                      ${hypre_mpi_flag}
    BUILD_COMMAND "$(MAKE)" all
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()
