# Truchas files in directory
#   drivers

set(ver_info ${Truchas_BINARY_DIR}/version_info.F90)
set(ver_info_in ${CMAKE_CURRENT_SOURCE_DIR}/drivers/version_info.F90.in)

add_custom_target(git_versioning ALL
  COMMAND ${CMAKE_COMMAND} "-DGIT_FOUND=${Git_FOUND}" "-DGIT=${GIT_EXECUTABLE}" "-DINFILE=${ver_info_in}" "-DOUTFILE=${ver_info}" -P "${CMAKE_SOURCE_DIR}/cmake/git_describe.cmake"
  BYPRODUCTS ${ver_info})


set(DRIVERS_SOURCE_FILES
           ${ver_info}
           drivers/drivers.F90
           drivers/physics_module.F90
           drivers/physical_constants.F90
           drivers/time_step_module.F90
           drivers/time_step_sync_type.F90
           drivers/simulation_event_queue.F90
           drivers/sim_event_queue_type.F90
           drivers/hijack_truchas.F90)


# Set compile flags
include(BuildWhitespaceString)
set(fc_flags
  -I${PGSLib_MODULE_DIR} -I${Truchas_utilities_dir})
build_whitespace_string(DRIVERS_COMPILE_FLAGS ${fc_flags})
set_source_files_properties(${DRIVERS_SOURCE_FILES} PROPERTIES
  COMPILE_FLAGS ${DRIVERS_COMPILE_FLAGS})

set_source_files_properties(drivers/drivers.F90 PROPERTIES
  COMPILE_DEFINITIONS "${Truchas_INFO_FLAGS}")

# Add the C source files
set(DRIVERS_C_SOURCE_FILES drivers/runinfo.c)
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DRIVERS_C_SOURCE_FILES})
list(APPEND Truchas_LIBRARY_SOURCE_FILES ${DRIVERS_SOURCE_FILES})
