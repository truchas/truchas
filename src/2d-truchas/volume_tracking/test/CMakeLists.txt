find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_executable(test_geometric_volume_tracker test_geometric_volume_tracker.F90)
target_link_libraries(test_geometric_volume_tracker volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(test_geometric_volume_tracker BEFORE PRIVATE
  "${VOLUME_TRACKING_2D_MODULE_DIR}")

add_executable(test_simple_volume_tracker test_simple_volume_tracker.F90)
target_link_libraries(test_simple_volume_tracker volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(test_simple_volume_tracker BEFORE PRIVATE
  "${VOLUME_TRACKING_2D_MODULE_DIR}")

add_test(NAME geometric_volume_tracker
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
    $<TARGET_FILE:test_geometric_volume_tracker> ${MPIEXEC_POSTFLAGS})
set_tests_properties(geometric_volume_tracker PROPERTIES
  PROCESSORS 1 LABELS "UNIT;2d-truchas")

add_test(NAME simple_volume_tracker
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
    $<TARGET_FILE:test_simple_volume_tracker> ${MPIEXEC_POSTFLAGS})
set_tests_properties(simple_volume_tracker PROPERTIES
  PROCESSORS 1 LABELS "UNIT;2d-truchas")

function(add_vof_test name script executable input_file gold_file)
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE} ${script}
      $<TARGET_FILE:${executable}> ${input_file} ${gold_file} ${MPIEXEC})
  set_tests_properties(${name} PROPERTIES
    PROCESSORS 4
    LABELS "INTEGRATED;2d-truchas"
    TIMEOUT 120)
endfunction()

add_subdirectory(advection)
add_subdirectory(vortex)
add_subdirectory(axisymmetric)
add_subdirectory(axial_flow)
add_subdirectory(generic)
add_subdirectory(truncated_volumes)
