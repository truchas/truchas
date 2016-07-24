if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR "ERROR: In-source builds are not allowed. Please "
    "create a separate directory and run cmake from there; for example,\n"
    "  $ mkdir MY_BUILD\n"
    "  $ cd MY_BUILD\n"
    "  $ cmake [OPTIONS] ..\n"
    "NB: You must first remove the CMakeCache.txt file and CMakeFiles "
    "directory that cmake just created before attempting to run cmake again; "
    "for example,\n"
    "  $ rm -r CMakeCache.txt CMakeFiles\n")
endif()
