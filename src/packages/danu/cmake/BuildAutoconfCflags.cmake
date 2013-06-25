include(PrintVariable)
include(SetMacros)
include(ListJoin)
include(CMakeParseArguments)

function(BUILD_AUTOCONF_CFLAGS OUTPUT)
  set(_oneValue "")
  set(_multiValue "FLAGS")
  set(_options "")
  cmake_parse_arguments(PARSE "${_options}" "${_oneValue}" "${_multiValue}" "${ARGN}")

  set(cflags ${CMAKE_C_FLAGS})
  if ( CMAKE_BUILD_TYPE )
    string(TOUPPER ${CMAKE_BUILD_TYPE} build_type) 
    append_set(cflags ${CMAKE_C_FLAGS_${build_type}})
  endif()

  append_set(cflags ${PARSE_FLAGS})
  list_join(VALUES ${cflags} GLUE " " OUTPUT _tmp_string)

  set(${OUTPUT} "${_tmp_string}" PARENT_SCOPE)


endfunction(BUILD_AUTOCONF_CFLAGS)
