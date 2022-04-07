if(SPHINX_FOUND)
  # Depends on importing the truchas python package at build,
  # which depends on truchas libraries being built.
  add_custom_target(reference-manual ALL
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${TruchasPython_BINARY_DIR}" ${SPHINX_EXECUTABLE} -b html
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/html
    DEPENDS TruchasPython gridmap fwrite
  )
endif()
