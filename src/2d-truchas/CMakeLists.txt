# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)
add_subdirectory(volume_tracking)

add_executable(mesh_example mesh_example.F90)

target_link_libraries(mesh_example truchas ${TRUCHAS_LINK_LIBS})
