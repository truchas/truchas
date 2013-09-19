# ############################################################################ #
#
# ADD_PYTRUCHAS_TEST(<test_name> <test_script>
#                    [CONFIGURATION [Debug|Release|...]]
#                    [DEPENDS test1 test2 ...] [LABELS lab1 lab2 ...]
#                    [ENVIRONMENT MYVAR1=VAL1 MYVAR2=VAR2 ... ]
#                    [PROCESSORS n] [RUN_SERIAL]
#                    [SERIAL_ONLY] [PARALLEL_ONLY]
#                    [TIMEOUT t] [WILL_FAIL])
#
# Add a TruchasTest script to the CTest suite. Test will be named <name> and
# <test_script> will be executed for the test. Other options are:
#
# CONFIGURATION [Debug|Release|...]       Only run this test under these 
#                                         configurations
# 
# DEPENDS test1 test2 ...                 Run test after test1, tes2,... complete.
#
# LABELS lab1 lab2 ...                    Labels associated with the test.
#
# ENVIRONMENT MYVAR1=VAL1 MYVAR2=VAL2 ... Set environment variables MYVAR1, 
#                                         MYVAR2 to VAL1, VAL2, ... during the
#                                         test. The PYTHONPATHS environment variable
#                                         is set in this macro to include the Truchas
#                                         Python build directory. Will append
#                                         additional paths if this is one of
#                                         the variables passed in.
#
# PROCESSORS n                            Number of processors required for
#                                         MPI (parallel) tests. Ignored if
#                                         ENABLE_MPI=False
#
# RUN_SERIAL                              Flag that forces CTest to NOT run other
#                                         tests in parallel while this test is
#                                         running.
#                                         
#
# SERIAL_ONLY                             Only run test if the Truchas binary
#                                         is a serial build
#
# PARALLEL_ONLY                           Only run test if the Truchas binary
#                                         is a MPI (parallel) build
#
# TIMEOUT t                               Test should be marked as timeout (FAIL)
#                                         if run time exceeds t seconds
#
# WILL_FAIL                               Test is expected to FAIL (!=0 exit)
#
# ############################################################################ #
include(CMakeParseArguments)

FUNCTION(ADD_PYTRUCHAS_TEST test_name test_script)

  # Parse the arguments
  set(options RUN_SERIAL SERIAL_ONLY PARALLEL_ONLY WILL_FAIL)
  set(oneValueArgs PROCESSORS TIMEOUT)
  set(multiValueArgs CONFIGURATION DEPENDS LABELS ENVIRONMENT)
  cmake_parse_arguments(MY_ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Exit now if PYTHON_EXECUTABLE has not been defined
  if ( NOT PYTHON_EXECUTABLE )
    message(FATAL_ERROR "Must define PYTHON_EXECUTABLE before calling ADD_PYTRUCHAS_TEST")
  endif()

  # Return immediately if this is a {SERAIL,PARALLEL}_ONLY and ENABLE_MPI
  # is not compatible.
  if ( ENABLE_MPI AND MY_ARG_SERIAL_ONLY )
    return()
  endif()

  if (NOT ENABLE_MPI AND MY_ARG_PARALLEL_ONLY)
    return()
  endif()

  # Now add the test
  if (IS_ABSOLUTE ${test_script})
    set(py_command ${PYTHON_EXECUTABLE} ${test_script} -v)
  else()
    set(py_command ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_script} -v)
  endif()  
  add_test(${test_name} ${py_command})

  # Set flag and one value properties
  set(test_properties)

  if ( MY_ARG_WILL_FAIL )
    list(APPEND test_properties WILL_FAIL TRUE)
  endif()

  if ( MY_ARG_PROCESSORS AND ENABLE_MPI )
    list(APPEND test_properties PROCESSORS ${MY_ARG_PROCESSORS})
  endif()

  if ( MY_ARG_RUN_SERIAL )
    list(APPEND test_properties RUN_SERIAL TRUE)
  endif()

  if ( MY_ARG_TIMEOUT )
    list(APPEND test_properties TIMEOUT ${MY_ARG_TIMEOUT})
  endif()

  if (test_properties)
    set_tests_properties(${test_name} PROPERTIES ${test_properties})
  endif()

  # Set the list properties, set_tests_properties does not work correctly 
  if ( MY_ARG_DEPENDS )
    set_property(TEST ${test_name} PROPERTY DEPENDS ${MY_ARG_DEPENDS})
  endif()

  if ( MY_ARG_LABELS )
    set_property(TEST ${test_name} PROPERTY LABELS ${MY_ARG_LABELS})
  endif()

  # Handle the PYTHONPATH environment variable

  # Treat this as a list then convert to a ':' string 
  set(py_paths)

  # Add Truchas Python and Danu Python build directories
  list(APPEND py_paths ${TruchasPython_BINARY_DIR} ${PyDanu_BINARY_DIR})


  # Preserve the user's PYTHONPATH
  if ( DEFINED ENV{PYTHONPATH} )
    string(REGEX REPLACE ":" ";" list_pythonpath "$ENV{PYTHONPATH}")
    list(APPEND py_paths ${list_pythonpath})
  endif()

  # Passed in as an argument? Remove it from the list
  if ( MY_ARG_ENVIRONMENT )
    set(idx 0)
    set(rm_idx)
    foreach ( var ${MY_ARG_ENVIRONMENT} )
      if ( "${var}" MATCHES "PYTHONPATH=" )
	string(REPLACE "PYTHONPATH=" "" var ${var})
	string(REPLACE ":" ";" var ${var})
	list(APPEND py_paths ${var})
	list(APPEND rm_idx ${idx})
      endif()
      math(EXPR idx "${idx}+1")
    endforeach()
    foreach(i ${rm_idx})
      list(REMOVE_AT MY_ARG_ENVIRONMENT ${i})
    endforeach()  
  endif()

  # Convert the list to a ":"
  if ( py_paths )
    list(GET py_paths 0 env_pypaths)
    list(REMOVE_AT py_paths 0)
    foreach(p ${py_paths})
      set(env_pypaths "${env_pypaths}:${p}")
    endforeach()
  endif()

  # Add to ENVIRONMENT arg and set the property
  if (MY_ARG_ENVIRONMENT)
    list(APPEND MY_ARG_ENVIRONMENT "PYTHONPATH=${env_pypaths}")
  else()
    set(MY_ARG_ENVIRONMENT "PYTHONPATH=${env_pypaths}")
  endif()  
  set_property(TEST ${test_name} PROPERTY ENVIRONMENT ${MY_ARG_ENVIRONMENT})

  get_test_property(${test_name} ENVIRONMENT testme)

ENDFUNCTION(ADD_PYTRUCHAS_TEST test_name test_script)

