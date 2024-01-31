# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

add_library(truchas-extra)
add_subdirectory(volume_fractions)
add_subdirectory(materials)
target_link_libraries(truchas-extra PUBLIC truchas)
target_include_directories(truchas-extra PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(volume_tracking)
add_subdirectory(heat_transfer)

add_executable(mesh_example mesh_example.F90)

target_link_libraries(mesh_example truchas ${TRUCHAS_LINK_LIBS})
