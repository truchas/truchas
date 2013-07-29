project(TBrookParser)

message(STATUS "configuring the tbrookparser directory")

# Add the Python and NumPy include directories
include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})

add_subdirectory(PythonPackages/TBrookParser/Extensions)

# Scripts to install
set(TB_SOURCE_SCRIPTS TBrookParse.py)

# Process the scripts
set(TB_SCRIPTS)
foreach(script ${TB_SOURCE_SCRIPTS})
  get_filename_component(base ${script} NAME)
  set(output_target ${TBrookParser_BINARY_DIR}/${base})
  configure_file(${script} ${output_target} @ONLY)
  list(APPEND TB_SCRIPTS ${output_target})
endforeach()

install(PROGRAMS ${TB_SCRIPTS}
        DESTINATION ${Truchas_BIN_INSTALL_DIR})
