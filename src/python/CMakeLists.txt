project(TruchasPython NONE)

set(PYTHON_FILES
  setup.py
  truchas/__init__.py
  truchas/TruchasData.py
  truchas/TruchasEnvironment.py
  truchas/TruchasTest.py
  fortran_write/__init__.py
  fortran_write/fortran_write.py
  bin/write-xdmf.py
  bin/write-gmv.py
  bin/write-restart.py
  bin/write-probe.py
  )

set(PYTHON_CONFIG_FILE truchas/TruchasConfig.py)

# preprocess the configuration file
get_target_property(Truchas_EXECUTABLE truchas_exe OUTPUT_NAME)
configure_file(${PYTHON_CONFIG_FILE}.in ${PYTHON_CONFIG_FILE} @ONLY)

# build & install rules
# copy .py files to the build directory
set(TruchasPython_BINARY_DIR ${CMAKE_CURRENT_BUILD_DIR})
foreach(f ${PYTHON_FILES})
  add_custom_command(OUTPUT ${f}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${f}
    DEPENDS ${f}
    )
endforeach(f)
add_custom_target(TruchasPython ALL DEPENDS ${PYTHON_FILES})

set(TruchasPython_INSTALL_PREFIX
  lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
install(FILES ${PYTHON_FILES} ${PYTHON_CONFIG_FILE}
  DESTINATION ${TruchasPython_INSTALL_PREFIX})

# build fortran write library
add_library(fwrite fortran_write/fortran_write.f90)
install(TARGETS fwrite
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
