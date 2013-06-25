function(ADD_PYTHON_TEST testname pyfile)

  set(_options "")
  set(_oneValueArgs "")
  set(_multiValueArgs "PYTHON_PATH")
  cmake_parse_arguments(PARSE "${_options}" "${_oneValueArgs}" "${_multiValueArgs}" ${ARGN} ) 

  # --- Locate the Python interpeter if PYTHON_EXECUTABLE is not set
  if ( NOT PYTHON_EXECUTABLE )
    find_package(PythonInterp)
    if ( NOT PYTHONINTERP_FOUND )
      message(FATAL_EROR "Attempted to add a Python Unit Test but failed to find Python")
    endif()  
  endif()  

  # --- Set the PYTHONPATH
  set(test_python_path "${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}")
  foreach(py_path ${PARSE_PYTHON_PATH})
    set(test_python_path "${py_path}:${test_python_path}")
  endforeach()
  
  # --- Build the command
  set(py_command  "${PYTHON_EXECUTABLE}" "${pyfile}" "-v")

  # --- Create the test
  add_test(${testname} ${py_command})

  # --- Set the environment for the test
  set_tests_properties(${testname} PROPERTIES
                       ENVIRONMENT "PYTHONPATH=${test_python_path}")

endfunction(ADD_PYTHON_TEST)
