# ############################################################################ #
#
#
#
# Find NumPy (Python extension) CMake Module
# ############################################################################ #
include(FindPackageHandleStandardArgs)

if ( NOT NUMPY_INCLUDE_DIRS )

  find_package(Python REQUIRED QUIET)

  
  # --- Determine the include path
  set(_tmp_py_file ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/findNumPyinclude.py)
  file(WRITE ${_tmp_py_file}
             "try: import numpy; print numpy.get_include()\nexcept: pass\n")
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${_tmp_py_file}"
                  RESULT_VARIABLE result
		  OUTPUT_VARIABLE output
		  ERROR_VARIABLE output
		  OUTPUT_STRIP_TRAILING_WHITESPACE
		  ERROR_STRIP_TRAILING_WHITESPACE)

  if ( result )
    message(SEND_ERROR "Failed to successfully import NumPy module")
    message(SEND_ERROR "Output from test script:\n${output}")

    set(NUMPY_INCLUDE_DIRS NUMPY_INCLUDE_DIRS-NOTFOUND)

  else()

    set(NUMPY_INCLUDE_DIRS ${output} ${PYTHON_INCLUDE_DIRS})

  endif()

endif()

if ( NOT NUMPY_VERSION )

  find_package(Python REQUIRED QUIET)

  # --- Now find the NumPy version
  set(_tmp_py_file ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/findNumPyversion.py)
  file(WRITE ${_tmp_py_file}
             "try: import numpy; print numpy.__version__\nexcept: pass\n")
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" "${_tmp_py_file}"
                  RESULT_VARIABLE result
		  OUTPUT_VARIABLE output
		  ERROR_VARIABLE output
		  OUTPUT_STRIP_TRAILING_WHITESPACE)

  if ( result )
    message(SEND_ERROR "Failed to successfully import NumPy module to find version")
    message(SEND_ERROR "Output from test script:\n${output}")

    set(NUMPY_VERSION NUMPY_VERSION-NOTFOUND)

  else()

    set(NUMPY_VERSION ${output})

  endif()

endif()

find_package_handle_standard_args(NumPy
                                  DEFAULT_MSG
				  NUMPY_VERSION
				  NUMPY_INCLUDE_DIRS)



