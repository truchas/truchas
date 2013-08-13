# ############################################################################ #
#
#  ADD_PYTHON_TEST(testname py_module 
#                  [PYTHON_PATH path1 [path2] [path3] ...]
#                 )
#
#  Add test testname with Python executable py_module. PYTHON_PATH is a list
#  of paths to add to the PYTHONPATH variable. By default the current binary
#  (build) directory is always added to the PYTHONPATH variable. The Python
#  script or module must run the appropriate test(s) when executed by
#  PYTHON_EXECUTABLE. Once the test is added, other properties such as LABEL
#  PROCESSORS can be added with set_test_properties.
#
# ############################################################################ #
function(ADD_PYTHON_TEST testname pyfile)

  set(_options "")
  set(_oneValueArgs "")
  set(_multiValueArgs "PYTHON_PATH")
  cmake_parse_arguments(PARSE "${_options}" "${_oneValueArgs}" "${_multiValueArgs}" ${ARGN} ) 

  # --- Fail here if PYTHON_EXECUTABLE is not defined
  if ( NOT PYTHON_EXECUTABLE )
    message(FATAL_ERROR "Python executable is not defined (PYTHON_EXECUTABLE). "
                        "Can not call add_python_test without this definition.")
  endif()		      

  # --- Set the PYTHONPATH , preserve the current PYTHONPATH
  if ( "$ENV{PYTHONPATH}" )
    set(test_python_path "${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}")
  else()  
    set(test_python_path ${CMAKE_CURRENT_BINARY_DIR})
  endif()  
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
