project(TruchasPython)

message(STATUS "configuring the Truchas Python directory")

# --- Python defintions and usefule macros
find_package(Python REQUIRED)
include(UsePython)

# --- Install

# Location 
if (NOT TruchasPython_INSTALL_PREFIX)
  set(version_short python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
  set(TruchasPython_INSTALL_PREFIX lib/${version_short}/site-packages)
endif()

# Variable (list) that holds files to install
set(TruchasPython_INSTALL_FILES)

# --- Define scripts
set(TruchasPython_COMPILE_SCRIPT ${Truchas_SOURCE_DIR}/scripts/PythonCompile.py)
set(TruchasPython_CONFIGURE_FILES ${TruchasPython_SOURCE_DIR}/cmake/configure-files.cmake)
set(TruchasPython_CONFIGURE_INCLUDE ${TruchasPython_SOURCE_DIR}/cmake/configure-include.cmake.in)

# --- Configure the include file

# Executable
get_target_property(Truchas_EXECUTABLE truchas_exe OUTPUT_NAME)

# Executable build location
print_variable(TruchasExe_BINARY_DIR)

# MPI flag
include(BoolEval)
bool_eval(ENABLE_MPI_FLAG ${ENABLE_MPI})

# Need the test directory, may not be added to the build yet.
if ( NOT TruchasTest_SOURCE_DIR)
  set(TruchasTest_SOURCE_DIR ${Truchas_SOURCE_DIR}/test)
endif()

if ( NOT TruchasTest_BINARY_DIR)
  set(TruchasTest_BINARY_DIR ${Truchas_BINARY_DIR}/test)
endif()


# configure
configure_file(${TruchasPython_CONFIGURE_INCLUDE}
               ${TruchasPython_BINARY_DIR}/configure-include.cmake
	       @ONLY)

# --- Add the directories
add_subdirectory(Truchas)

# --- Copy Python files
set(TruchasPython_FILES
                 FileFilter
		 TruchasTest)

# Each target is <file> with out the extension
foreach(f ${TruchasPython_FILES})
  set(mod_target ${f})
  set(in_file ${TruchasPython_SOURCE_DIR}/${f}.py.in)
  set(py_file ${TruchasPython_BINARY_DIR}/${f}.py)
  set(pyc_file ${py_file}c)
  set(config_args -DINFILE:STRING=${in_file} -DOUTFILE:STRING=${py_file})
  add_custom_command(OUTPUT ${py_file} ${pyc_file}
                     DEPENDS ${in_file}
                     COMMAND
		     ${CMAKE_COMMAND} ${config_args} -P ${TruchasPython_CONFIGURE_FILES}
		     COMMAND
		     ${PYTHON_EXECUTABLE} ${TruchasPython_COMPILE_SCRIPT} ${py_file}
		     COMMENT "Processing ${f}")
  add_custom_target(${mod_target} ALL DEPENDS ${py_file} ${pyc_file})	
  list(APPEND TruchasPython_INSTALL_FILES ${py_file} ${pyc_file})
endforeach()		  
                 
install(PROGRAMS ${TruchasPython_INSTALL_FILES}
        DESTINATION ${TruchasPython_INSTALL_PREFIX})

install(PROGRAMS sample-truchas.py.in
        DESTINATION bin)