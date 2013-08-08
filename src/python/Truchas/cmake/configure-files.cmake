# ############################################################################ #
# 
# CMake script to configure the Python files
#
# Usage:
#  cmake -DINFILE:STRING=<in_file> -DOUTFILE:STRING=<outfile> -P configure-files.cmake
# ############################################################################ #


# Include file contains important variables
set(INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/configure-include.cmake)
if ( EXISTS ${INCLUDE_FILE} )
  include(${INCLUDE_FILE})
else()
  message(FATAL_ERROR "Include file ${INCLUDE_FILE} does not exist")
endif()

if(NOT INFILE)
  message(FATAL_ERROR "Input file is not defined")
endif()

if(NOT OUTFILE)
  message(FATAL_ERROR "Output file is not defined")
endif()

configure_file(${INFILE} ${OUTFILE} @ONLY)


