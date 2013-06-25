# - Find SWIG
# This module finds an installed SWIG.  It sets the following variables:
#  SWIG_FOUND - set to true if SWIG is found
#  SWIG_DIR - the directory where swig is installed
#  SWIG_EXECUTABLE - the path to the swig executable
#  SWIG_VERSION   - the version number of the swig executable
#
# The minimum required version of SWIG can be specified using the
# standard syntax, e.g. FIND_PACKAGE(SWIG 1.1)
#
# All information is collected from the SWIG_EXECUTABLE so the
# version to be found can be changed from the command line by
# means of setting SWIG_EXECUTABLE
#

#=============================================================================
# Copyright 2004-2009 Kitware, Inc.
# Copyright 2011 Mathieu Malaterre <mathieu.malaterre@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# FULL CMake Copyright Notice
#CMake - Cross Platform Makefile Generator
#Copyright 2000-2011 Kitware, Inc., Insight Software Consortium
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#
#* Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
#* Neither the names of Kitware, Inc., the Insight Software Consortium,
#  nor the names of their contributors may be used to endorse or promote
#  products derived from this software without specific prior written
#  permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#------------------------------------------------------------------------------
#
#The above copyright and license notice applies to distributions of
#CMake in source and binary form.  Some source files contain additional
#notices of original copyright by their contributors; see each source
#for details.  Third-party software packages supplied with CMake under
#compatible licenses provide their own copyright notices documented in
#corresponding subdirectories.
#
#------------------------------------------------------------------------------
#
#CMake was initially developed by Kitware with the following sponsorship:
#
# * National Library of Medicine at the National Institutes of Health
#   as part of the Insight Segmentation and Registration Toolkit (ITK).
#
# * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#   Visualization Initiative.
#
# * National Alliance for Medical Image Computing (NAMIC) is funded by the
#   National Institutes of Health through the NIH Roadmap for Medical Research,
#   Grant U54 EB005149.
#
# * Kitware, Inc.
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

FIND_PROGRAM(SWIG_EXECUTABLE NAMES swig2.0 swig)

IF(SWIG_EXECUTABLE)

    IF((NOT SWIG_VERSION) AND (NOT SWIG_DIR) )
    EXECUTE_PROCESS(COMMAND ${SWIG_EXECUTABLE} -swiglib
      OUTPUT_VARIABLE SWIG_swiglib_output
      ERROR_VARIABLE SWIG_swiglib_error
      RESULT_VARIABLE SWIG_swiglib_result)

    IF(SWIG_swiglib_result) 
      IF(SWIG_FIND_REQUIRED)
        MESSAGE(SEND_ERROR "Command \"${SWIG_EXECUTABLE} -swiglib\" failed with output:\n${SWIG_swiglib_error}")
      ELSE(SWIG_FIND_REQUIRED)
        MESSAGE(STATUS "Command \"${SWIG_EXECUTABLE} -swiglib\" failed with output:\n${SWIG_swiglib_error}")
      ENDIF(SWIG_FIND_REQUIRED)
    ELSE(SWIG_swiglib_result)
      STRING(REGEX REPLACE "[\n\r]+" ";" SWIG_swiglib_output ${SWIG_swiglib_output})
      # force the path to be computed each time in case SWIG_EXECUTABLE has changed.
      SET(SWIG_DIR SWIG_DIR-NOTFOUND)
      FIND_PATH(SWIG_DIR swig.swg PATHS ${SWIG_swiglib_output})
      IF(SWIG_DIR)
        SET(SWIG_USE_FILE ${CMAKE_ROOT}/Modules/UseSWIG.cmake)
        EXECUTE_PROCESS(COMMAND ${SWIG_EXECUTABLE} -version
          OUTPUT_VARIABLE SWIG_version_output
          ERROR_VARIABLE SWIG_version_output
          RESULT_VARIABLE SWIG_version_result)
        IF(SWIG_version_result)
          MESSAGE(SEND_ERROR "Command \"${SWIG_EXECUTABLE} -version\" failed with output:\n${SWIG_version_output}")
        ELSE(SWIG_version_result)
          STRING(REGEX REPLACE ".*SWIG Version[^0-9.]*\([0-9.]+\).*" "\\1"
            SWIG_version_output "${SWIG_version_output}")
          SET(SWIG_VERSION ${SWIG_version_output} CACHE STRING "Swig version" FORCE)
        ENDIF(SWIG_version_result)
      ENDIF(SWIG_DIR)
    ENDIF(SWIG_swiglib_result)
  ENDIF((NOT SWIG_VERSION) AND (NOT SWIG_DIR))  
ENDIF(SWIG_EXECUTABLE)

#INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SWIG  REQUIRED_VARS SWIG_EXECUTABLE SWIG_DIR
                                        VERSION_VAR SWIG_VERSION )
