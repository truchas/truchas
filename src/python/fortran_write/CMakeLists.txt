project(FortranWrite NONE)

set(PYTHON_SCRIPT_FILES
  __init__.py
  fortran_write.py
  )

# build rules
foreach(f ${PYTHON_SCRIPT_FILES})
  add_custom_command(OUTPUT ${f}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${f}
    DEPENDS ${f}
    )
endforeach(f)
add_custom_target(FortranWritePython ALL DEPENDS ${PYTHON_SCRIPT_FILES})

add_library(fwrite fortran_write.f90)

# install rules
foreach(f ${PYTHON_SCRIPT_FILES})
  install(FILES ${FortranWrite_BINARY_DIR}/${f}
    DESTINATION ${TruchasPython_INSTALL_PREFIX}/fortran_write)
endforeach(f)

install(TARGETS fwrite
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
