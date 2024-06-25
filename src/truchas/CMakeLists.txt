project(TruchasExe C Fortran)

# hack for bad clang optimization on apple
if(APPLE)
  add_compile_definitions(TRUCHAS_ALLOW_UNSAFE_VECTORIZATION)
endif()

# Add the SOURCE_DIR and BINARY_DIR to the include directories
set(CMAKE_INCLUDE_CURRENT_DIR True)

# --- Library

set(ver_info ${Truchas_BINARY_DIR}/version_info.F90)
set(ver_info_in ${CMAKE_CURRENT_SOURCE_DIR}/version_info.F90.in)

add_custom_target(git_versioning ALL
  COMMAND ${CMAKE_COMMAND} "-DGIT_FOUND=${Git_FOUND}" "-DGIT=${GIT_EXECUTABLE}" "-DINFILE=${ver_info_in}" "-DOUTFILE=${ver_info}" -P "${CMAKE_SOURCE_DIR}/cmake/git_describe.cmake"
  BYPRODUCTS ${ver_info})


# Cannot set the compile definitions for these individual files from the
# CMakeLists.txt file in the respective subdirectories for some reason.

set_property(SOURCE drivers/drivers.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "${Truchas_INFO_FLAGS}")

if(CMAKE_SYSTEM_NAME MATCHES Linux)
  set_property(SOURCE utilities/get_process_size.c APPEND PROPERTY
    COMPILE_DEFINITIONS LINUX)
endif()

add_library(truchas ${ver_info})
set_target_properties(truchas PROPERTIES LINKER_LANGUAGE Fortran)

# Library sources (and unit tests)
add_subdirectory(communication)
add_subdirectory(discrete_operators)
add_subdirectory(distributed_mesh)
add_subdirectory(drivers)
add_subdirectory(exodus)
add_subdirectory(functions)
add_subdirectory(mesh_functions)
add_subdirectory(data_mapping)
add_subdirectory(input)
add_subdirectory(materials)
add_subdirectory(ustruc)
add_subdirectory(ode)
add_subdirectory(output)
add_subdirectory(partitioning)
add_subdirectory(flow)
add_subdirectory(solid_mechanics)
add_subdirectory(physics)
add_subdirectory(setup)
add_subdirectory(solver)
add_subdirectory(utilities)
add_subdirectory(toolpath)
add_subdirectory(truchasio)
add_subdirectory(vtkhdf)

target_link_libraries(truchas PUBLIC
    LAPACK::LAPACK exodus hypre petaca hdf5::hdf5 metis scorpio MPI::MPI_Fortran netcdf
    ${CMAKE_DL_LIBS} # to pick up the dlopen and dlclose symbols
)
if(USE_PORTAGE)
  target_link_libraries(truchas PUBLIC portage::portage)
endif()

target_include_directories(truchas PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utilities>
    $<BUILD_INTERFACE:${Truchas_MODULE_DIR}>
)

# --- Main executable

add_executable(truchas_exe drivers/main.F90)
target_link_libraries(truchas_exe truchas)

# Keep binary name consistent with old naming convention
set(truchas_exe_name t-${Truchas_BUILD_ID})
set_target_properties(truchas_exe PROPERTIES
                      OUTPUT_NAME ${truchas_exe_name}
                      EchoString "Built Truchas executable ${truchas_exe_name}")

# --- Install rules
install(TARGETS truchas truchas_exe
        RUNTIME DESTINATION ${Truchas_BIN_INSTALL_DIR}
        LIBRARY DESTINATION ${Truchas_LIBRARY_INSTALL_DIR}
        ARCHIVE DESTINATION ${Truchas_LIBRARY_INSTALL_DIR})
