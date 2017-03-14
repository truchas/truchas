set(SRC
    scorpio_c_binding.F90
    scorpio_file_type.F90
    truchas_h5_outfile.F90
    scorpio_ext.c
    )

add_library(truchas_io ${SRC})
target_link_libraries(truchas_io PUBLIC truchas::danu truchas::pgslib_fe
  truchas::pgslib_mpi truchas::scorpio)

add_library(truchas::io ALIAS truchas_io)
set_property(TARGET truchas_io APPEND PROPERTY INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(TARGETS truchas_io
        EXPORT  truchas_io
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

