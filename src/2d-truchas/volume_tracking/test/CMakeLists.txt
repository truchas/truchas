find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(VOF_TEST_SUPPORT_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")
add_library(vof_2d_test_support OBJECT
  support/read_inputfile.F90
  support/vof_2d_test_driver.F90
  support/gaussian_quadrature_vofinit.F90
)
set_target_properties(vof_2d_test_support PROPERTIES
  Fortran_MODULE_DIRECTORY "${VOF_TEST_SUPPORT_MODULE_DIR}")
target_link_libraries(vof_2d_test_support volume_tracking_2d)
# The 2D tracker modules intentionally share names with mainline's 3D
# tracker modules. Compile the test support against the 2D modules.
target_include_directories(vof_2d_test_support BEFORE PRIVATE
  "${VOLUME_TRACKING_2D_MODULE_DIR}")

function(add_vof_test_executable target source)
  add_executable(${target} ${source} $<TARGET_OBJECTS:vof_2d_test_support>)
  target_link_libraries(${target} volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
  target_include_directories(${target} BEFORE PRIVATE
    "${VOF_TEST_SUPPORT_MODULE_DIR}"
    "${VOLUME_TRACKING_2D_MODULE_DIR}")
endfunction()

function(add_vof_legacy_test name executable input_file reference_file)
  add_test(NAME ${name}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:${executable}> ${input_file} ${reference_file} ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${name} PROPERTIES PROCESSORS 4 LABELS 2d-truchas)
endfunction()

function(add_vof_vtkhdf_test name script executable input_file reference_file)
  add_test(NAME ${name}_vtkhdf
    COMMAND ${Python3_EXECUTABLE} ${script}
      $<TARGET_FILE:${executable}> ${input_file} ${reference_file} ${MPIEXEC})
  set_tests_properties(${name}_vtkhdf PROPERTIES
    PROCESSORS 4
    LABELS "INTEGRATED;2d-truchas"
    TIMEOUT 120)
endfunction()

add_subdirectory(advection)
add_subdirectory(vortex)
add_subdirectory(axisymmetric)
add_subdirectory(axial_flow)
add_subdirectory(truncated_volumes)
