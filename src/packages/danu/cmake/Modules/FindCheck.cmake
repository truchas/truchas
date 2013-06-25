# -*- mode: cmake -*-

#
# Amanzi CHECK Find Module
#
# Usage:
#    Control the search through CHECK_DIR or setting environment variable
#    CHECK_ROOT to the CHECK installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CHECK_FOUND            (BOOL)       Flag indicating if CHECK was found
#    CHECK_INCLUDE_DIR      (PATH)       Path to the CHECK include file
#    CHECK_INCLUDE_DIRS     (LIST)       List of all required include files
#    CHECK_LIBRARY_DIR      (PATH)       Path to the CHECK library
#    CHECK_LIBRARY          (FILE)       CHECK library
#    CHECK_LIBRARIES        (LIST)       List of all required CHECK libraries
#
#    Additional variables
#    CHECK_VERSION          (STRING)     CHECK Version string
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

include(PrintVariable)


if ( CHECK_LIBRARIES AND CHECK_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(CHECK_LIBRARIES AND CHECK_INCLUDE_DIRS)

    # Cache variables
    if(CHECK_DIR)
        set(CHECK_DIR "${CHECK_DIR}" CACHE PATH "Path to search for CHECK include and library files")
    endif()

    if(CHECK_INCLUDE_DIR)
        set(CHECK_INCLUDE_DIR "${CHECK_INCLUDE_DIR}" CACHE PATH "Path to search for CHECK include files")
    endif()

    if(CHECK_LIBRARY_DIR)
        set(CHECK_LIBRARY_DIR "${CHECK_LIBRARY_DIR}" CACHE PATH "Path to search for CHECK library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) CHECK_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) CHECK_DIR/<include,include/CHECK++,include/check++>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(check_inc_names "check.h")
    if (CHECK_INCLUDE_DIR)

        if (EXISTS "${CHECK_INCLUDE_DIR}")

            find_path(test_include_path
                      NAMES ${check_inc_names}
                      HINTS ${CHECK_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT test_include_path)
                message(WARNING "Can not locate ${check_inc_names} in ${CHECK_INCLUDE_DIR}")
            endif()
            set(CHECK_INCLUDE_DIR "${test_include_path}")

        else()
            message(WARNING "CHECK_INCLUDE_DIR=${CHECK_INCLUDE_DIR} does not exist")
            set(CHECK_INCLUDE_DIR "CHECK_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(check_inc_suffixes "include")
        if(CHECK_DIR)

            if (EXISTS "${CHECK_DIR}" )

                find_path(CHECK_INCLUDE_DIR
                          NAMES ${check_inc_names}
                          HINTS ${CHECK_DIR}
                          PATH_SUFFIXES ${check_inc_suffixes}
                          NO_DEFAULT_PATH)
                 print_variable(CHECK_INCLUDE_PATH)
            else()
                 message(WARNING "CHECK_DIR=${CHECK_DIR} does not exist")
                 set(CHECK_INCLUDE_DIR "CHECK_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(CHECK_INCLUDE_DIR
                      NAMES ${check_inc_names}
                      PATH_SUFFIXES ${check_inc_suffixes})

        endif()

    endif()

    if ( NOT CHECK_INCLUDE_DIR )
        message(WARNING "Can not locate Check include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) CHECK_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) CHECK_DIR/<lib,lib/CHECK++,lib/check++>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(check_lib_names "check")
    if (CHECK_LIBRARY_DIR)

        if (EXISTS "${CHECK_LIBRARY_DIR}")

            find_library(CHECK_LIBRARY
                         NAMES ${check_lib_names}
                         HINTS ${CHECK_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(WARNING "CHECK_LIBRARY_DIR=${CHECK_LIBRARY_DIR} does not exist")
            set(CHECK_LIBRARY "CHECK_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND check_lib_suffixes "lib" "Lib")
        if(CHECK_DIR)

            if (EXISTS "${CHECK_DIR}" )

                find_library(CHECK_LIBRARY
                             NAMES ${check_lib_names}
                             HINTS ${CHECK_DIR}
                             PATH_SUFFIXES ${check_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(WARNING "CHECK_DIR=${CHECK_DIR} does not exist")
                 set(CHECKLIBRARY "CHECK_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(CHECK_LIBRARY
                         NAMES ${check_lib_names}
                         PATH_SUFFIXES ${check_lib_suffixes})

        endif()

    endif()

    if ( NOT CHECK_LIBRARY )
        message(WARNING "Can not locate Check library")
    endif()    

   
    # CHECK does not have any prerequisite libraries
    set(CHECK_INCLUDE_DIRS ${CHECK_INCLUDE_DIR})
    set(CHECK_LIBRARIES    ${CHECK_LIBRARY})

   
endif(CHECK_LIBRARIES AND CHECK_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(CHECK DEFAULT_MSG
                                           CHECK_LIBRARIES
                                           CHECK_INCLUDE_DIRS)

# Define the version
if ( CHECK_INCLUDE_DIR )
    set(check_h "${CHECK_INCLUDE_DIR}/check.h")
    file(STRINGS "${check_h}" check_version_string REGEX "^#define CHECK_VERSION ")
    string(REGEX REPLACE "^#define CHECK_VERSION ([0-9]+).*$" "\\1" check_version "${check_version_string}")
    set(CHECK_VERSION "${check_version}")
endif()    

mark_as_advanced(
  CHECK_VERSION
  CHECK_INCLUDE_DIR
  CHECK_INCLUDE_DIRS
  CHECK_LIBRARY
  CHECK_LIBRARIES
  CHECK_LIBRARY_DIR
)
