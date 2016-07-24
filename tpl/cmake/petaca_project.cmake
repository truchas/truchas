if(SEARCH_FOR_PETACA)
  message(STATUS "Searching for a suitable Petaca library ...")
  find_package(PETACA)
endif()

if(PETACA_FOUND)
  list(APPEND projects_found "Petaca")
else()
  list(APPEND projects_to_build "Petaca")
  get_filename_component(yajl_library_dir ${YAJL_LIBRARY} DIRECTORY)
  externalproject_add(petaca
    DEPENDS yajl
    PREFIX petaca
    URL ${TARFILE_DIR}/petaca-c49f1f1.tar.gz
    URL_MD5 d9b3a06674405635ae75045f05032594
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_Fortran_COMPILER:PATH=${CMAKE_Fortran_COMPILER}
               -D CMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
               -D CMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
               -D CMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
               -D YAJL_INCLUDE_DIR:PATH=${YAJL_INCLUDE_DIR}
               -D YAJL_LIBRARY_DIR:PATH=${yajl_library_dir}
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
               -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()
