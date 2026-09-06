find_package(Python3 COMPONENTS Interpreter REQUIRED)

function(add_vof_test name script input_file gold_file)
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE} ${script}
      $<TARGET_FILE:truchas_2d> ${input_file} ${gold_file} ${MPIEXEC})
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
