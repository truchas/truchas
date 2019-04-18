project(TruchasPython NONE)

# ZJJ 03-29-19
# We need to copy scripts from truchas/ and fortran_write/
# as-is. The truchas/TruchasConfig.py.in and bin/* files need
# to be configured with a copy and replace on @keywords@.
# unfortunately, CMake's configure_file happens at configure
# time rather than build time. The result is that changing
# these files causes not merely a copy-and-configure to the
# build directory, but a complete reconfigure and rebuild
# of the entire project. The workaround, if one wanted to
# use it, would be to put configure_file in a separate
# CMake script and passing it the needed variables through
# add_custom_command. See the mailing list discussion here:
#
# https://cmake.org/pipermail/cmake/2015-April/060364.html

set(PYTHON_SCRIPT_FILES
  truchas/__init__.py
  truchas/TruchasData.py
  truchas/TruchasEnvironment.py
  truchas/TruchasTest.py
  )

set(PYTHON_BIN_FILES
  bin/write-xdmf.py
  bin/write-gmv.py
  bin/write-restart.py
  bin/write-probe.py
  )

# Both Build and Install config files get configured
# and placed in the build tree, but only the Install
# config file is moved to the install directory.
set(PYTHON_CONFIG_FILES
  truchas/TruchasConfigBuild.py
  truchas/TruchasConfigInstall.py
  )

set(TruchasPython_BINARY_DIR ${CMAKE_CURRENT_BUILD_DIR})
set(TruchasPython_INSTALL_PREFIX
  ${Truchas_LIBRARY_INSTALL_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

# preprocess the bin and configuration files
get_target_property(Truchas_EXECUTABLE truchas_exe OUTPUT_NAME)
foreach(f ${PYTHON_BIN_FILES} ${PYTHON_CONFIG_FILES})
  configure_file(${f}.in ${f} @ONLY)
endforeach(f)

# build rules: copy .py files to the build directory
foreach(f ${PYTHON_SCRIPT_FILES})
  add_custom_command(OUTPUT ${f}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${f}
    DEPENDS ${f}
    )
endforeach(f)
add_custom_target(TruchasPython ALL DEPENDS ${PYTHON_SCRIPT_FILES})

foreach(f ${PYTHON_SCRIPT_FILES} truchas/TruchasConfigInstall.py)
  install(FILES ${TruchasPython_BINARY_DIR}/${f}
    DESTINATION ${TruchasPython_INSTALL_PREFIX}/truchas)
endforeach(f)

foreach(f ${PYTHON_BIN_FILES})
  install(PROGRAMS ${TruchasPython_BINARY_DIR}/${f}
    DESTINATION bin)
endforeach(f)


add_subdirectory(fortran_write)
