# ############################################################################ #
#
#
#   DANU CMake 
#      Python Interface Source File Directory
#
# ############################################################################ #
if (ENABLE_Python)
  project(PyDanu C)

# ---------------------------------------------------------------------------- #
# Python Interfaces built with SWIG www.swig.org
# ---------------------------------------------------------------------------- #
  if(NOT TARGET ${SWIG_BUILD_TARGET}) # SWIG is part of the build, no search
    find_package(SWIG)
    if ( NOT SWIG_FOUND )
      message(FATAL_ERROR "Can not locate SWIG. Can not build Python interfaces")
    endif()  
  endif()  

  include(UseSWIG)

# ---------------------------------------------------------------------------- #
# Include directories 
# ---------------------------------------------------------------------------- #
  # Add External Package Includes
  include_directories(${HDF5_INCLUDE_DIRS})
  include_directories(${PYTHON_INCLUDE_DIRS})
  include_directories(${NUMPY_INCLUDE_DIRS})

  # Danu Includes
  include_directories("${Danu_SOURCE_DIR}/src")
  include_directories(${PyDanu_SOURCE_DIR})

# ---------------------------------------------------------------------------- #
# Compiler flags
# ---------------------------------------------------------------------------- #
  set(EXTRA_LINK_FLAGS)

  if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
      set(EXTRA_LINK_FLAGS "-undefined dynamic_lookup")
    endif()
  endif()

# ---------------------------------------------------------------------------- #
# SWIG flags
# ---------------------------------------------------------------------------- #

  set(SWIG_FLAGS "-Wall -v")
  set_property(GLOBAL PROPERTY CMAKE_SWIG_FLAGS ${SWIG_FLAGS})


# ---------------------------------------------------------------------------- #
# SWIG file definitions
# ---------------------------------------------------------------------------- #
  # SWIG Interface Files
  set(SWIG_I_FILES numpy.i
                   Danu.i
                   H5Obj.i
                   File.i
		   Output.i
		   Mesh.i
		   Simulations.i
                   exception_handling.i
                   hid_t.i)

  # C code for the interfaces
  set(SWIG_C_FILES except.c
                   PyUtils.c
                   NumPyUtils.c
                   H5Object.c
                   File.c
		   Output.c
		   Mesh.c
		   Simulations.c)

  # SWIG flags passed through a property
  set_source_files_properties(${SWIG_I_FILES}
                              PROPERTIES 
                              CMAKE_SWIG_FLAGS "${SWIG_FLAGS}")             

  # Create the Python module Danu
  # Danu.i is the interface file that defines the
  # Python library
  SWIG_ADD_MODULE(Danu python Danu.i ${SWIG_C_FILES})
  SWIG_LINK_LIBRARIES(Danu danu)
  set_target_properties(_Danu
                        PROPERTIES
                        LINK_FLAGS "${EXTRA_LINK_FLAGS}")
  if ( TARGET ${SWIG_BUILD_TARGET} )
    add_dependencies(_Danu ${SWIG_BUILD_TARGET})
  endif()

# ---------------------------------------------------------------------------- #
# GMV 
# ---------------------------------------------------------------------------- #

# Create a Python library that interfaces with the GMV library
# General Mesh Viewer (GMV)
set(GMV_C_SOURCE GMV.c gmvwrite.c except.c )
set(GMV_SWIG_I_FILES numpy.i exception_handling.i)

set_source_files_properties(${GMV_SIWG_I_FILES}
                              PROPERTIES 
                              CMAKE_SWIG_FLAGS "${SWIG_FLAGS}")

# The gmvwrite.c file includes malloc.h
# On some systems (MacOSX) this has been deprecated
# and malloc is now part of stdlib.h
# Check to see if the the include file is present.
# If not make sure malloc exists and exit if it is
# not found.
include(CheckIncludeFile)
check_include_file(malloc.h HAVE_MALLOC_H)
if ( NOT HAVE_MALLOC_H )
  include(CheckFunctionExists)
  check_function_exists(malloc HAVE_MALLOC)
  if ( NOT HAVE_MALLOC )
    message(FATAL_ERROR "C function malloc is not found. gmvwrite.c requires malloc")
  endif()
else()
  set_source_files_properties(gmvwrite.c PROPERTIES
                              COMPILE_DEFINITIONS "HAVE_MALLOC_H")
endif()			   


SWIG_ADD_MODULE(GMV python GMV.i ${GMV_C_SOURCE})
set_target_properties(_GMV
                      PROPERTIES
                      LINK_FLAGS "${EXTRA_LINK_FLAGS}")
if (TARGET ${SWIG_BUILD_TARGET} )		    
  add_dependencies(_GMV ${SWIG_BUILD_TARGET})
endif()

# ---------------------------------------------------------------------------- #
# Danu Python install path
# ---------------------------------------------------------------------------- #

  # Modules to install
  set(PyDanu_MODULES
            ${PyDanu_BINARY_DIR}/Danu.py
            ${PyDanu_BINARY_DIR}/_Danu.so
            ${PyDanu_BINARY_DIR}/GMV.py
            ${PyDanu_BINARY_DIR}/_GMV.so
            ${PyDanu_SOURCE_DIR}/xdmf_conversion.py)

  install(FILES ${PyDanu_MODULES} 
          DESTINATION ${Danu_Python_INSTALL_DIR})

  # Scripts to install
  set(PyDanu_SOURCE_SCRIPTS
      ${PyDanu_SOURCE_DIR}/danu_report.py
      ${PyDanu_SOURCE_DIR}/truchas-gmv-parser.py
      ${PyDanu_SOURCE_DIR}/xdmf-parser.py)

    # Process the scripts  
  set(PyDanu_SCRIPTS)
  foreach( script ${PyDanu_SOURCE_SCRIPTS} )
    get_filename_component(base ${script} NAME)
    set(output_target ${CMAKE_CURRENT_BINARY_DIR}/${base})
    configure_file(${script} ${output_target} @ONLY)
    list(APPEND PyDanu_SCRIPTS ${output_target})
  endforeach()  

  install(PROGRAMS ${PyDanu_SCRIPTS} 
          DESTINATION ${Danu_BIN_INSTALL_DIR})



# ---------------------------------------------------------------------------- #
# Danu Python Unit Tests
# ---------------------------------------------------------------------------- #

#DEBUG  if ( ENABLE_Tests )
#DEBUG
#DEBUG    enable_testing()
#DEBUG
#DEBUG    include(AddPythonTest)
#DEBUG
#DEBUG    # Removing the leading directories here makes the REGEX easier
#DEBUG    file(GLOB PyDanu_TEST_MODULES RELATIVE ${PyDanu_SOURCE_DIR}/test "${PyDanu_SOURCE_DIR}/test/*.py")
#DEBUG
#DEBUG    foreach ( py_module ${PyDanu_TEST_MODULES} )
#DEBUG      get_filename_component(py_module_we ${py_module} NAME_WE)
#DEBUG      string(REGEX REPLACE "Test" "py" test_name "${py_module_we}")
#DEBUG      add_python_test(${test_name} ${PyDanu_SOURCE_DIR}/test/${py_module})
#DEBUG
#DEBUG    endforeach()  
#DEBUG
#DEBUG  endif()  
endif()    


  


