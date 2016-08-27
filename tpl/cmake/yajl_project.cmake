# Builds both shared and static library; not controllable.

if(SEARCH_FOR_YAJL)
  message(STATUS "Searching for a suitable YAJL library ...")
  find_package(YAJL "2.0.4")
endif()

if(YAJL_FOUND)
  list(APPEND projects_found "YAJL")
else()
  list(APPEND projects_to_build "YAJL")
  set(YAJL_VERSION "2.1.0")
  externalproject_add(yajl
    PREFIX yajl
    URL ${TARFILE_DIR}/yajl-${YAJL_VERSION}.tar.gz
    URL_MD5 6887e0ed7479d2549761a4d284d3ecb0
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
               -D CMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
  set(YAJL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
  if(BUILD_SHARED_LIBS)
    set(YAJL_LIBRARY "${CMAKE_INSTALL_PREFIX}/libyajl${CMAKE_SHARED_LIBRARY_SUFFIX}")
  else()
    set(YAJL_LIBRARY "${CMAKE_INSTALL_PREFIX}/libyajl_s${CMAKE_STATIC_LIBRARY_SUFFIX}")
  endif()
  set(YAJL_INCLUDE_DIRS "${YAJL_INCLUDE_DIR}")
  set(YAJL_LIBRARIES "${YAJL_LIBRARY}")
endif()
