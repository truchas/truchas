project(TBrookParser)

message(STATUS "configuring the tbrookparser directory")

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

install(FILES ${TB_SCRIPTS}
        DESTINATION ${Truchas_BIN_INSTALL_DIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
