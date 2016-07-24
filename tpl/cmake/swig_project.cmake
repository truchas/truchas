if(SEARCH_FOR_SWIG)
  message(STATUS "Searching for a suitable swig executable ...")
  find_package(SWIG)
  if(SWIG_FOUND)
    if(SWIG_VERSION VERSION_LESS "2" OR NOT (SWIG_VERSION VERSION_LESS "3"))
      set(SWIG_FOUND FALSE)
      message(STATUS "Could NOT find a suitable SWIG: version 2 is required")
    endif()
  endif()
endif()

if(SWIG_FOUND)
  list(APPEND projects_found "SWIG")
else()
  list(APPEND projects_to_build "SWIG")
  set(SWIG_VERSION "2.0.4")
  set(SWIG_EXECUTABLE "${CMAKE_INSTALL_PREFIX}/bin/swig")
  externalproject_add(swig
    PREFIX swig
    URL ${TARFILE_DIR}/swig-${SWIG_VERSION}.tar.gz
    URL_MD5 4319c503ee3a13d2a53be9d828c3adc0
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
                      CC=${CMAKE_C_COMPILER}
                      CXX=${CMAKE_CXX_COMPILER}
                      CFLAGS=${cflags}
                      CPPFLAGS=${cppflags}
                      --without-pcre
                      --with-python=${PYTHON_EXECUTABLE}
                      --prefix=${CMAKE_INSTALL_PREFIX}
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

