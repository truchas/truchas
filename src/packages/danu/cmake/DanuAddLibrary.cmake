# ########################################################################### #
#
# Danu Add Library Function
#
# ########################################################################### #


function(DANU_ADD_LIBRARY library_name)

  include(CMakeParseArguments)
  include(PrintVariable)

  set(options   "STATIC_ONLY;NO_INSTALL;NO_INSTALL_HEADERS")
  set(oneValue  "DESTINATION;HEADER_DESTINATION")
  set(arguments "SOURCE;HEADERS;LINK_LIBS;DEPENDENCIES")
  cmake_parse_arguments(ARG "${options}" "${oneValue}" "${arguments}" "${ARGN}")

  list(LENGTH ARG_HEADERS dummy)

  if(NOT library_name)
    message(FATAL_ERROR "Must define a library target name")
  endif()

  if(NOT ARG_SOURCE)
    message(FATAL_ERROR "Must define a source file list for ${library_name}")
  endif()

  # Create the target for the library
  if(ARG_STATIC_ONLY)
    add_library(${library_name} STATIC ${ARG_SOURCE} ${ARG_HEADERS})
  else()  
    add_library(${library_name} ${ARG_SOURCE} ${ARG_HEADERS})
  endif()  

  # Add link libraries
  if(ARG_LINK_LIBS)
    target_link_libraries(${library_name} ${ARG_LINK_LIBS})
  endif()

  # Add dependencies
  if(ARG_DEPENDENCIES)
    add_dependencies("${library_name}" "${ARG_DEPENDENCIES}")
  endif()

  # Create the install target
  if ( NOT ARG_NO_INSTALL )
    set(lib_dest "lib")
    if(ARG_DESTINATION)
      set(lib_dest ${ARG_DESTINATION})
    endif()
    install(TARGETS ${library_name} EXPORT ${library_name} DESTINATION ${lib_dest})  
  endif()  

  # Add headers files to the install target
  if(ARG_HEADERS AND (NOT ARG_NO_INSTALL_HEADERS))
    set(inc_dest include)
    if ( ARG_HEADER_DESTINATION )
      set(inc_dest ${ARG_HEADER_DESTINATION})
    endif()  
    install(FILES ${ARG_HEADERS} DESTINATION ${inc_dest})
  endif( ARG_HEADERS AND (NOT ARG_NO_INSTALL_HEADERS))  
    

endfunction(DANU_ADD_LIBRARY)
