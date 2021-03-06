include(CMakeParseArguments)
#function(LIST_JOIN VALUES GLUE OUTPUT)
function(LIST_JOIN)
  set(_oneValue "GLUE;OUTPUT")
  set(_multiValue "VALUES")
  set(_options "")
  cmake_parse_arguments(PARSE "${_options}" "${_oneValue}" "${_multiValue}" "${ARGN}")
  string (REGEX REPLACE "([^\\]|^);" "\\1${PARSE_GLUE}" _TMP_STR "${PARSE_VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${PARSE_OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

