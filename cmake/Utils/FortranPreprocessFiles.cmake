include(CMakeParseArguments)
# - FORTRAN_PREPROCESS_FILES
# 
# Usage:
#
# FORTRAN_PREPROCESS_FILES(<output variable>
#                          FILES file1 [file2 [file3 [...]]]
#                          FPP_EXECUTABLE <preprocessor>
#                          [FPP_FLAGS flag1 flag2 ... ]
#                          [GREP_FILTERS arg1 arg2 arg3]
#                          [PROCESS_TARGET <target>] 
#)
#                                                     
MACRO(_define_filenames infile outfile fullname)
  get_filename_component(base "${infile}" NAME_WE)  
  get_filename_component(ext  "${infile}" EXT)
  get_filename_component(abspath "${infile}" ABSOLUTE)

  set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}.f90)
  set(fullname ${CMAKE_CURRENT_SOURCE_DIR}/${infile})

ENDMACRO()

FUNCTION(FORTRAN_PREPROCESS_FILES mylist)

  # Parse the arguments
  set(options)
  set(oneValueArgs PROCESS_TARGET)
  set(multiValueArgs FILES FPP_EXECUTABLE GREP_FILTERS FPP_FLAGS)
  cmake_parse_arguments(MY "${options}" "${oneValueArgs}" "${multiValueArgs}" 
    ${ARGN})

  # Check the output name
  if (NOT MY_FILES)
    message(FATAL_ERROR "Must provide a list of files to process")
  endif()

  # Set the optional values
  if (NOT MY_FPP_EXECUTABLE )
    message(FATAL_ERROR "Must provide a preprocessor executable")
  endif()

  # Build the grep filter pipes
  # Always grep out the lines with # at the beginning
  set(grep_filter grep -v ^\#)
  foreach(a ${MY_GREP_FILTERS})
    list(APPEND grep_filter | grep -v "${a}")
  endforeach()

  # Loop through each file and add a custom command
  set(outfiles)
  foreach(infile ${MY_FILES})
    _define_filenames(${infile} outfile fullname)
    #message(STATUS "outfile=${outfile}")
    #message(STATUS "fullname=${fullname}")
    add_custom_command(OUTPUT "${outfile}"
        COMMAND ${MY_FPP_EXECUTABLE} ${MY_FPP_FLAGS} ${fullname} | ${grep_filter} > ${outfile}
	DEPENDS ${fullname}
	COMMENT "Preprocessing ${infile}"
	VERBATIM)
    list(APPEND outfiles ${outfile})  
  endforeach()


  # Now create a target for this command
  if (MY_PROCESS_TARGET)
    add_custom_target(${MY_PROCESS_TARGET} ALL DEPENDS ${outfiles})
  endif()

  # Now push the output file list  up to the parent 
  set("${mylist}" ${outfiles} PARENT_SCOPE)



ENDFUNCTION()
