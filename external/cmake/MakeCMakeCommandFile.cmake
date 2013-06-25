# #  -*- mode: cmake -*-

#
# MAKE_CMAKE_COMMAND_FILE(cmake_file
#                         SCRIPT script_file
#                         ACTION string
#                         EP_NAME project_name
#                         WORK_DIRECTORY dir)
#
include(CMakeParseArguments)
FUNCTION(MAKE_CMAKE_COMMAND_FILE cmake_file)

  # Parse arguments
  set(_flags    "")
  set(_oneValue SCRIPT EP_NAME WORK_DIRECTORY ACTION)
  set(_multiValue "")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})

  set(_script   ${PARSE_SCRIPT})
  set(_ep_name  ${PARSE_EP_NAME})
  set(_work_dir ${PARSE_WORK_DIRECTORY})
  set(_action   ${PARSE_ACTION})

  file(TO_CMAKE_PATH ${CMAKE_CURRENT_LIST_DIR}/../templates template_dir)
  set(template_file  ${template_dir}/general-command-step.cmake.in)

  configure_file(${template_file} ${cmake_file} @ONLY@)


ENDFUNCTION()  
