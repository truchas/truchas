if(SPHINX_FOUND)
  add_custom_target(reference-manual ALL
      ${SPHINX_EXECUTABLE} -b html
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}/html
  )
endif()
