# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

if(CMAKE_BUILD_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" TRUCHAS_2D_BUILD_TYPE_UPPER)
  set(TRUCHAS_2D_BUILD_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${TRUCHAS_2D_BUILD_TYPE_UPPER}}")
else()
  set(TRUCHAS_2D_BUILD_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()
string(STRIP "${TRUCHAS_2D_BUILD_FLAGS}" TRUCHAS_2D_BUILD_FLAGS)
configure_file(simulation_build_info.F90.in simulation_build_info.F90 @ONLY)

add_library(truchas-extra
  ${CMAKE_CURRENT_BINARY_DIR}/simulation_build_info.F90
  simulation_log_type.F90
  simulation_environment_type.F90
  simulation_class.F90
  simulation_command_line_type.F90
  simulation_provenance.F90)
add_subdirectory(mesh)
add_subdirectory(material_distribution)
add_subdirectory(materials)
target_link_libraries(truchas-extra PUBLIC truchas)
target_include_directories(truchas-extra PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(Test)
add_subdirectory(physics/volume_tracking)
add_subdirectory(physics/flow)
add_subdirectory(physics/thermal)
add_subdirectory(simulations/flow)
add_subdirectory(simulations/thermal)
add_subdirectory(simulations/flow_thermal)
add_subdirectory(simulations/vof)

# Existing test trees will be reorganized separately.
add_subdirectory(volume_tracking/test)
add_subdirectory(flow/Test)
add_subdirectory(heat_transfer/test)
add_subdirectory(nonisothermal/test)

add_executable(truchas_2d truchas_2d_main.F90)
set_target_properties(truchas_2d PROPERTIES OUTPUT_NAME truchas-2d)
target_link_libraries(truchas_2d nsht2d flow_sim2d thermal_sim2d vof2d truchas-extra truchas ${TRUCHAS_LINK_LIBS})
