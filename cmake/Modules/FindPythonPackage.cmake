# Check if a python package is available, send fatal error if it isn't,
# or if the found version is less than a given minimum.

function(find_python_package PACKAGE MIN_VERSION)
  if(NOT DEFINED PYTHON_EXECUTABLE)
    message(STATUS "${PYTHON_EXECUTABLE}")
    message(FATAL_ERROR "Python executable not defined. Could not search for ${PACKAGE}")
  endif()

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
    "import ${PACKAGE}; print(${PACKAGE}.__version__, end='')"
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE found_package_version
    ERROR_QUIET)

  if(exit_code)
    message(FATAL_ERROR "Python module ${PACKAGE} not found")
  elseif(found_package_version VERSION_LESS MIN_VERSION)
    message(FATAL_ERROR "Could NOT find Python module ${PACKAGE}: Found \
unsuitable version \"${found_package_version}\", but required is at least \"${MIN_VERSION}\"")
  else()
    message(STATUS "Found Python module ${PACKAGE} (found version ${found_package_version})")
  endif()
endfunction(find_python_package)
