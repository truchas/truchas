# ############################################################################ #
#
#
#   DANU CMake 
#    
#      Create an imported target
#
# ############################################################################ #
include(ParseArguments)

include(PrintVariable)

function(CREATE_IMPORT_TARGET)

    # Macro: _print_usage
    macro(_print_usage)
        message("CREATE_IMPORT_TARGET(file_name>\n"
                "         TARGET_NAME target_name\n"
                "         [STATIC_ONLY] (ignored if file_name is executable\n"
                "         [EXECUTABLE]  (file_name is executable not library\n")
    endmacro(_print_usage)     

    # Read in args
    set(options "EXECUTABLE;STATIC_ONLY")
    set(arg_name "TARGET_NAME")
    parse_arguments(IMPORT "${arg_name}" "${options}" "${ARGN}")
    set(IMPORT_FILE "${IMPORT_DEFAULT_ARGS}")
    if ( NOT IMPORT_FILE )
        _print_usage()
        message(FATAL_ERROR "Must define file name")
    endif()
    if ( NOT IMPORT_TARGET_NAME )
        _print_usage()
        message(FATAL_ERROR "Failed to define a target name")
    endif()
    if (TARGET "${IMPORT_TARGET_NAME}")
        message(FATAL_ERROR "Can not create an import target ${IMPORT_TARGET_NAME}."
                            " It already exists.")
    endif()

    get_filename_component(IMPORT_FILE_FULLPATH ${IMPORT_FILE} ABSOLUTE)


    if (IMPORT_EXECUTABLE)
        add_executable("${IMPORT_TARGET_NAME}" IMPORTED)
    else()
        if(IMPORT_STATIC_ONLY)
            add_library("${IMPORT_TARGET_NAME}" STATIC IMPORTED)
        else()
            add_library("${IMPORT_TARGET_NAME}" UNKNOWN IMPORTED)
        endif()
    endif()     
    
    set_target_properties("${IMPORT_TARGET_NAME}"
                          PROPERTIES
                          IMPORTED_LOCATION "${IMPORT_FILE_FULLPATH}")


endfunction(CREATE_IMPORT_TARGET)             
