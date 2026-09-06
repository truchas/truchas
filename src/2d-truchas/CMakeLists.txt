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
configure_file(src/core/simulation_build_info.F90.in simulation_build_info.F90 @ONLY)

add_library(truchas-extra
  ${CMAKE_CURRENT_BINARY_DIR}/simulation_build_info.F90
  src/core/simulation_log_type.F90
  src/core/simulation_environment_type.F90
  src/core/simulation_class.F90
  src/core/simulation_command_line_type.F90
  src/core/simulation_provenance.F90)
add_subdirectory(src/mesh)
add_subdirectory(src/material_distribution)
add_subdirectory(src/materials)
target_link_libraries(truchas-extra PUBLIC truchas)
target_include_directories(truchas-extra PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(src/core/test)
add_subdirectory(src/physics/volume_tracking)
add_subdirectory(src/physics/flow)
add_subdirectory(src/physics/thermal)
add_subdirectory(src/simulations/flow)
add_subdirectory(src/simulations/thermal)
add_subdirectory(src/simulations/flow_thermal)
add_subdirectory(src/simulations/vof)

add_subdirectory(src/physics/volume_tracking/test)
add_subdirectory(src/physics/flow/test)
add_subdirectory(src/physics/thermal/test)
add_subdirectory(src/simulations/flow/test)
add_subdirectory(src/simulations/thermal/test)
add_subdirectory(src/simulations/flow_thermal/test)
add_subdirectory(src/simulations/vof/test)

add_executable(truchas_2d src/truchas_2d_main.F90)
set_target_properties(truchas_2d PROPERTIES OUTPUT_NAME truchas-2d)
target_link_libraries(truchas_2d nsht2d flow_sim2d thermal_sim2d vof2d truchas-extra truchas ${TRUCHAS_LINK_LIBS})
