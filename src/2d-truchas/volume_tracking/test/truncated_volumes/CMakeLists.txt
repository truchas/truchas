add_executable(trunc_vol truncvolume_2dtest_driver.F90)
target_link_libraries(trunc_vol volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(trunc_vol BEFORE PRIVATE "${VOLUME_TRACKING_2D_MODULE_DIR}")

add_test(NAME truncated_volumes
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:trunc_vol> ${MPIEXEC_POSTFLAGS})
set_tests_properties(truncated_volumes PROPERTIES PROCESSORS 1 LABELS 2d-truchas)
