# Check if a python package is available, send fatal error if it isn't.
# Should add handling for minimum versions and REQUIRED flag.

function(find_python_package PACKAGE)
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
  else()
    message(STATUS "Found Python module ${PACKAGE} (found version ${found_package_version})")
  endif()
endfunction(find_python_package)
