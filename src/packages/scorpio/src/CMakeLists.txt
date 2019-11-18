cmake_minimum_required(VERSION 3.0)

project(scorpio)

# MPI dependency
add_library(p::mpi INTERFACE IMPORTED)
set_property(TARGET p::mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    ${MPI_C_INCLUDE_PATH})
set_property(TARGET p::mpi PROPERTY INTERFACE_LINK_LIBRARIES
    ${MPI_C_LIBRARIES})

# HDF5 dependency
add_library(p::hdf5 INTERFACE IMPORTED)
set_property(TARGET p::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    ${HDF5_INCLUDE_DIRS})
set_property(TARGET p::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES
    ${HDF5_C_LIBRARIES})

set(SRC
    scorpio.c
    scorpioftn.c
    scorpio_read.c
    scorpio_util.c
    scorpio_write.c
    )

add_library(scorpio ${SRC})
target_compile_definitions(scorpio PRIVATE -DLINUX_FTN)
target_link_libraries(scorpio PUBLIC p::hdf5 p::mpi)
add_library(truchas::scorpio ALIAS scorpio)

set_property(TARGET scorpio APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(TARGETS scorpio
        EXPORT  scorpio
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(example)
