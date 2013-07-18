# UPDATE_COMPILER_FLAGS(LANGUAGE FLAGS flag1 [flag2 ...] [PREAPPEND] [IGNORE_REPEAT])
#
# Update the CMAKE_<LANGUAGE>_FLAGS variable
#   LANGUAGE              Language (C, CXX or Fortran)
#   FLAGS                 List of flags to add
#   PREAPPEND             Flag that indicates flags are preappended to the variable
#   IGNORE_REPEAT         By default, macro will not add a flag already present
#                         this flag overrides that behavior
#
MACRO(UPDATE_COMPILER_FLAGS lang)

  # parse the arguments
  include(CMakeParseArguments)
  set(options PREAPPEND IGNORE_REPEAT)
  set(oneValueArgs "")
  set(multiValueArgs FLAGS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Zero out the _oldflags and _newflags. This is a 
  # macro these variables are defined in the calling space.
  set(_oldflags "")
  set(_newflags "")

  # CMAKE_<lang>_FLAGS is a string, want a list 
  if (CMAKE_${lang}_FLAGS)
    if(UNIX)
      separate_arguments(_oldflags UNIX_COMMAND "${CMAKE_${lang}_FLAGS}")
    elseif(WIN32)
      separate_arguments(_oldflags WINDOWS_COMMAND "${CMAKE_${lang}_FLAGS}")
    endif()
  endif()

  # Now add the new flags
  if (ARG_PREAPPEND)
    list(APPEND _newflags ${ARG_FLAGS})
    list(APPEND _newflags ${_oldflags})
  else()
    foreach(_flag ${_oldflags})
      string(STRIP ${_flag} _strip)
      list(APPEND _newflags ${_strip})
    endforeach()
    list(APPEND _newflags ${ARG_FLAGS})
  endif()

  # get rid of repeated flags
  if (NOT ARG_IGNORE_REPEATS)
    list(REMOVE_DUPLICATES _newflags)
  endif()

  # Update the CMAKE_<lang>_FLAGS variable
  list(GET _newflags 0 _first)
  set(CMAKE_${lang}_FLAGS  ${_first})
  list(REMOVE_AT _newflags 0)
  foreach(_flag ${_newflags})
    set(CMAKE_${lang}_FLAGS "${CMAKE_${lang}_FLAGS} ${_flag}")
  endforeach()
  
ENDMACRO(UPDATE_COMPILER_FLAGS lang)
