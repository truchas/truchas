project(GridMapping NONE)

set(PYTHON_SCRIPT_FILES
  __init__.py
  grid_mapping.py
  )

# build rules
foreach(f ${PYTHON_SCRIPT_FILES})
  add_custom_command(OUTPUT ${f}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${f}
    DEPENDS ${f}
    )
endforeach(f)
add_custom_target(GridMappingPython ALL DEPENDS ${PYTHON_SCRIPT_FILES})

if(NOT DOCUMENTATION_ONLY)
  include_directories(${Truchas_MODULE_DIR})
  add_library(gridmap grid_mapping_wrapper.F90)
  add_dependencies(gridmap truchas ${Truchas_DEPENDS})
  target_link_libraries(gridmap truchas)
endif()

# install rules
foreach(f ${PYTHON_SCRIPT_FILES})
  install(FILES ${GridMapping_BINARY_DIR}/${f}
    DESTINATION ${TruchasPython_INSTALL_PREFIX}/grid_mapping)
endforeach(f)

if(NOT DOCUMENTATION_ONLY)
  install(TARGETS gridmap
    LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
    ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
endif()
