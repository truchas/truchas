set(SRC
    truchasio.F90
    io_scorpio.c
    )

add_library(truchas_io ${SRC})
target_link_libraries(truchas_io PUBLIC truchas::danu truchas::pgslib_fe
  truchas::pgslib_mpi truchas::scorpio)

add_library(truchas::io ALIAS truchas_io)
set_property(TARGET truchas_io APPEND PROPERTY INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
