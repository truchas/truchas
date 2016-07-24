project(TruchasHarness NONE)

# --- Variable (list) to hold all the install files
set(TruchasHarness_INSTALL_FILES)

# --- Python files sans extension
set(TruchasHarness_FILES
                  __init__
		  __main__
		  main
		  binary
		  error
		  mpi
		  output)

# --- Setup the configure include file
configure_file(${TruchasPython_CONFIGURE_INCLUDE}
               ${TruchasHarness_BINARY_DIR}/configure-include.cmake
	       @ONLY)

# --- Process the files
foreach(f ${TruchasHarness_FILES})
  set(in_file ${TruchasHarness_SOURCE_DIR}/${f}.py.in)
  set(py_file ${TruchasHarness_BINARY_DIR}/${f}.py)
  set(pyc_file ${py_file}c)
  set(config_args -DINFILE:STRING=${in_file} -DOUTFILE:STRING=${py_file})
  add_custom_command(OUTPUT ${py_file} ${pyc_file}
                     DEPENDS ${in_file}
                     COMMAND
		     ${CMAKE_COMMAND} ${config_args} -P ${TruchasPython_CONFIGURE_FILES}
		     COMMAND
		     ${PYTHON_EXECUTABLE} ${TruchasPython_COMPILE_SCRIPT} ${py_file}
		     COMMENT "Processing ${f}")
  add_custom_target(Truchas${f} ALL DEPENDS ${py_file} ${pyc_file})		   
  list(APPEND TruchasHarness_INSTALL_FILES ${py_file} ${pyc_file})
endforeach()  

# --- Install rules
install(PROGRAMS ${TruchasHarness_INSTALL_FILES}
         DESTINATION ${TruchasPython_INSTALL_PREFIX}/Truchas)
                 
