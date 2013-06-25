# - Find python interpreter

# CMake has split finding the Python executable and libraries. 
# This FindPython module combines the functionality of both and
# sets the version string.
# This code sets the following variables:
# 
#  PYTHONINTERP_FOUND        - Was the Python executable found
#  PYTHON_EXECUTABLE         - path to the Python interpreter
#  PYTHON_VERSION_STRING     - Version of the executable
#  PYTHONLIBS_FOUND          - have the Python libs been found
#  PYTHON_LIBRARIES          - path to the python library
#  PYTHON_include_PATH       - path to where Python.h is found (deprecated)
#  PYTHON_include_DIRS       - path to where Python.h is found
#  PYTHON_DEBUG_LIBRARIES    - path to the debug library
#

#=============================================================================
# Copyright 2005-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# FULL CMake Copyright notice and license information
#Copyright 2000-2009 Kitware, Inc., Insight Software Consortium. All rights
#reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#
#Redistributions of source code must retain the above copyright notice, this 
#list of conditions and the following disclaimer.
#
#Redistributions in binary form must reproduce the above copyright notice, 
#this list of conditions and the following disclaimer in the documentation 
#and/or other materials provided with the distribution.
#
#Neither the names of Kitware, Inc., the Insight Software Consortium, nor 
#the names of their contributors may be used to endorse or promote products 
#derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#POSSIBILITY OF SUCH DAMAGE.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)
include(FindPackageHandleStandardArgs)
include(PrintVariable)


# Find the python executable...lifted from CMake FindPythonInterp module
FIND_PROGRAM(PYTHON_EXECUTABLE
  NAMES python2.7 python2.6 python2.5 python2.4 python2.3 python2.2 python2.1 python2.0 python1.6 python1.5 python
  PATHS
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.5\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.4\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.3\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.2\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.1\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.0\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.6\\InstallPath]
  [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\1.5\\InstallPath]
  )

if(NOT PYTHON_EXECUTABLE)
  if(Python_FIND_REQUIRED)
    message(FATAL_ERROR "Can not locate Python executable")
  endif()
endif()

# Determine the Python version
execute_process(COMMAND
                ${PYTHON_EXECUTABLE} "-V"
                OUTPUT_VARIABLE python_version_out
                ERROR_VARIABLE  python_version_out
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)
string(REPLACE "Python" "" python_version_str ${python_version_out}) 
string(STRIP ${python_version_str} PYTHON_VERSION_STRING)


# Find the major, minor and patch versions
if ( PYTHON_VERSION_STRING )
  set(PYTHON_VERSION ${PYTHON_VERSION_STRING})
  string(REPLACE "." ";" _tmp_list "${PYTHON_VERSION_STRING}")
  list(LENGTH _tmp_list _tmp_len)
  if ( _tmp_list AND  ( "${_tmp_len}" GREATER "1" ) )
    list(GET _tmp_list 0 PYTHON_VERSION_MAJOR)
    list(GET _tmp_list 1 PYTHON_VERSION_MINOR)
    list(GET _tmp_list 2 PYTHON_VERSION_PATCH)
  endif()
else()
  print_variable(python_version_out)
endif()


# Once the executable has been found, we now search for the
# include and library directories using this binary.
# Use the executable to determine the include and libraries
# paths. No INDENTING! Messses up python command.
# Lifted from http://code.google.com/p/numexpr/source/browse/FindPythonLibsNew.cmake
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
   "from distutils import sysconfig as s;import sys;import struct;
print(s.PREFIX);
print(s.get_python_inc(plat_specific=True));
print(s.get_python_lib(plat_specific=True));
"
   RESULT_VARIABLE _python_success
   OUTPUT_VARIABLE _python_out
   ERROR_VARIABLE  _python_err
   OUTPUT_STRIP_TRAILING_WHITESPACE)

# Blow up here if the above failed
if (NOT _python_success MATCHES 0 )
  if (Python_FIND_REQUIRED)
    message(FATAL_ERROR "Can not locate Python installation information: ${_python_err}")
  endif()
endif()

# Turn output into a list, easier to access
string(REGEX REPLACE "[\n\r]+" ";" _python_out ${_python_out})

# Define the include directory
list(GET _python_out 1 _inc_search_dir)
find_path(PYTHON_include_DIR
          Python.h
          PATHS ${_inc_search_dir}
          NO_DEFAULT_PATH)

# Define the Python library
list(GET _python_out 0 PYTHON_PREFIX)
find_library(PYTHON_LIBRARY
             NAMES "python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
             PATHS ${PYTHON_PREFIX}
	     PATH_SUFFIXES lib lib64 libs
	     NO_DEFAULT_PATH)

# For backwards compatibility
set(PYTHON_INCLUDE_DIRS "${PYTHON_include_DIR}")
set(PYTHON_LIBRARIES    "${PYTHON_LIBRARY}")

# handle the QUIETLY and REQUIRED arguments and set PYTHON_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Python
                                  DEFAULT_MSG
                                  PYTHON_EXECUTABLE
                                  PYTHON_VERSION
                                  PYTHON_LIBRARIES
                                  PYTHON_INCLUDE_DIRS)

mark_as_advanced(PYTHON_EXECUTABLE)
mark_as_advanced(PYTHON_VERSION_STRING)
mark_as_advanced(PYTHON_VERSION_MAJOR)
mark_as_advanced(PYTHON_VERSION_MINOR)
mark_as_advanced(PYTHON_VERSION_PATCH)
mark_as_advanced(PYTHON_include_DIRS)
mark_as_advanced(PYTHON_LIBRARIES)
