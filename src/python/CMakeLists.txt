project(TruchasPython)

message(STATUS "configuring the Truchas Python directory")

# --- Python defintions and usefule macros
find_package(Python REQUIRED)
include(UsePython)

# --- Install

# Location
if (NOT TruchasPython_INSTALL_PREFIX)
  set(version_short python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
  set(TruchasPython_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib/${version_short}/site-packages)
endif()

# Variable (list) that holds files to install
set(TruchasPython_INSTALL_FILES)

# --- Define scripts
set(TruchasPython_COMPILE_SCRIPT ${Truchas_SOURCE_DIR}/scripts/PythonCompile.py)

# --- Add the directories
add_subdirectory(Truchas)

# --- Copy Python files
set(TruchasPython_MODULES
                 FileFilter.py
		 TruchasTest.py)

# Each target is <file> with out the extension
foreach(module ${TruchasPython_MODULES})
  get_filename_component(module_we ${module} NAME_WE)
  set(mod_target ${module_we})
  set(old_file ${TruchasPython_SOURCE_DIR}/${module})
  set(new_file ${TruchasPython_BINARY_DIR}/${module})
  add_custom_command(OUTPUT ${new_file} ${new_file}c
		    COMMENT "Generating ${module}"
		    DEPENDS ${old_file}
		    COMMAND
                    ${CMAKE_COMMAND} -E copy_if_different ${old_file} ${new_file}
		    COMMAND
                    ${PYTHON_EXECUTABLE} ${TruchasPython_COMPILE_SCRIPT} ${new_file})
  add_custom_target(${mod_target} ALL DEPENDS ${new_file} ${new_file}c)		  
  list(APPEND TruchasPython_INSTALL_FILES ${new_file} ${new_file}c)
endforeach()		  
                 
install(PROGRAMS ${TruchasPython_INSTALL_FILES}
        DESTINATION ${TruchasPython_INSTALL_PREFIX})

install(PROGRAMS sample-truchas.py.in
        DESTINATION bin)
