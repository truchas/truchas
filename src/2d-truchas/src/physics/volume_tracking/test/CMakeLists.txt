include(CTest)

foreach(target IN ITEMS test_t2d_geometric_volume_tracker test_t2d_simple_volume_tracker)
  add_executable(${target} ${target}.F90)
  target_link_libraries(${target} volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
  add_test(NAME ${target}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${target} PROPERTIES
    PROCESSORS 1 LABELS "UNIT;VOLUME_TRACKING;2d-truchas")
endforeach()

add_subdirectory(truncated_volumes)
