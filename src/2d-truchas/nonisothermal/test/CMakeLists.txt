include(CTest)

add_executable(test_ns_ht_2d_enthalpy_advector test_ns_ht_2d_enthalpy_advector.F90)
target_link_libraries(test_ns_ht_2d_enthalpy_advector nsht2d truchas-extra truchas ${Truchas_LINK_LIBS})

add_executable(test_ns_ht_2d_solver test_ns_ht_2d_solver.F90)
target_link_libraries(test_ns_ht_2d_solver nsht2d truchas-extra truchas ${Truchas_LINK_LIBS})

foreach(nproc IN ITEMS 1 4)
  add_test(NAME ns_ht_2d_enthalpy_advector-${nproc}p
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
      $<TARGET_FILE:test_ns_ht_2d_enthalpy_advector> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ns_ht_2d_enthalpy_advector-${nproc}p PROPERTIES
    PROCESSORS ${nproc}
      LABELS "UNIT;NONISOTHERMAL;2d-truchas")

  add_test(NAME ns_ht_2d_solver-${nproc}p
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
      $<TARGET_FILE:test_ns_ht_2d_solver> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ns_ht_2d_solver-${nproc}p PROPERTIES
    PROCESSORS ${nproc}
    LABELS "UNIT;NONISOTHERMAL;2d-truchas")
endforeach()

find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_python_package(vtk 9.0)
add_test(NAME ns_ht_2d_channel
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/channel/test.py
    $<TARGET_FILE:ns_ht_2d_main>
    ${CMAKE_CURRENT_SOURCE_DIR}/channel/input.json ${MPIEXEC} 0.0)
set_tests_properties(ns_ht_2d_channel PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)

add_test(NAME ns_ht_2d_channel_rotated
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/channel/test.py
    $<TARGET_FILE:ns_ht_2d_main>
    ${CMAKE_CURRENT_SOURCE_DIR}/channel/rotated-input.json ${MPIEXEC} 30.0)
set_tests_properties(ns_ht_2d_channel_rotated PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)

add_test(NAME ns_ht_2d_nat_conv
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/nat-conv/test.py
    $<TARGET_FILE:ns_ht_2d_main> ${MPIEXEC})
set_tests_properties(ns_ht_2d_nat_conv PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 600)

add_test(NAME ns_ht_2d_hydrostatic
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/hydrostatic/test.py
    $<TARGET_FILE:ns_ht_2d_main> ${MPIEXEC})
set_tests_properties(ns_ht_2d_hydrostatic PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)

add_test(NAME ns_ht_2d_inviscid_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/inviscid/test.py
          $<TARGET_FILE:ns_ht_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/inviscid/input.json
          ${MPIEXEC})
set_tests_properties(ns_ht_2d_inviscid_channel PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)

add_test(NAME ns_ht_2d_inviscid_channel_rotated
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/inviscid/test.py
          $<TARGET_FILE:ns_ht_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/inviscid/rotated-input.json
          ${MPIEXEC} 30.0)
set_tests_properties(ns_ht_2d_inviscid_channel_rotated PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)

add_test(NAME ns_ht_2d_multifluid_channel
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/multifluid_channel/test.py
          $<TARGET_FILE:ns_ht_2d_main>
          ${CMAKE_CURRENT_SOURCE_DIR}/multifluid_channel/input.json
          ${MPIEXEC})
set_tests_properties(ns_ht_2d_multifluid_channel PROPERTIES
  PROCESSORS 4
  LABELS "INTEGRATION;NONISOTHERMAL;2d-truchas"
  TIMEOUT 120)
