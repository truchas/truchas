

include(PrintVariable)

macro(APPEND_SET VARIABLE)
  set(${VARIABLE} ${${VARIABLE}} ${ARGN})
endmacro(APPEND_SET)

macro(APPEND_PATH FILE_LIST)

  include(CMakeParseArguments)
  set(options   "")
  set(oneValue  "PATH")
  set(arguments "FILES")
  cmake_parse_arguments(ARG "${options}" "${oneValue}" "${arguments}" "${ARGN}")

  foreach( _file ${ARG_FILES} )
    set(_tmp ${ARG_PATH}/${_file})
    append_set(${FILE_LIST} ${_tmp})
  endforeach()  

endmacro(APPEND_PATH)
  

macro(DISTRIBUTION_DEFINITIONS)

  message(STATUS "Defining package targets")

  # Version information
  set(CPACK_PACKAGE_VERSION_MAJOR "${Danu_MAJOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_MINOR "${Danu_MINOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_PATCH "${Danu_PATCH_LEVEL}")
  set(CPACK_PACKAGE_VERSION       "${Danu_VERSION}")

  # Basic package information
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Danu provides basic interfaces to an HDF5 simulation output file format")
  set(CPACK_PACKAGE_VENDOR "Los Alamos National Laboratory")

  #set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)

  # Tar file definitions
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "danu-${Danu_VERSION}" CACHE INTERNAL "Tarfile basename")

  # Files to ignore 
  set(CPACK_SOURCE_IGNORE_FILES
      "^${Danu_SOURCE_DIR}/.hg"
      "^${Danu_SOURCE_DIR}/build.*"
      "^${Danu_SOURCE_DIR}/examples"
      "^${Danu_SOURCE_DIR}/test/test_*"
      ".*.swp"
      "Makefile"
      "CMakeFiles"
      "CMakeCache.txt"
      "cmake_install.cmake"
      )

  # This adds the 'package' target to the build system
  include(CPack)

endmacro(DISTRIBUTION_DEFINITIONS)

macro(ADD_CHECK_TEST test_src_file)

  get_filename_component(file_we ${test_src_file} NAME_WE)
  set(test_executable "${file_we}.x")
  string(REGEX REPLACE "check_" "" test_name "${file_we}")
  #print_variable(test_executable)
  #print_variable(test_name)


  add_executable(${test_executable} ${test_src_file})
  target_link_libraries(${test_executable} ${ARGN})
  if ( BUILD_Check )
    add_dependencies(${test_executable} ${Check_BUILD_TARGET})
  endif()
  add_test(${test_name} ${test_executable})

endmacro(ADD_CHECK_TEST)

macro(ADD_FRUIT_TEST test_src_file)

  get_filename_component(file_we ${test_src_file} NAME_WE)
  set(test_executable "${file_we}.x")
  string(REGEX REPLACE "funit_" "fort_" test_name "${file_we}")
  #print_variable(test_executable)
  #print_variable(test_name)


  add_executable(${test_executable} ${test_src_file})
  target_link_libraries(${test_executable} ${ARGN})
  add_test(${test_name} ${test_executable})

endmacro(ADD_FRUIT_TEST)
