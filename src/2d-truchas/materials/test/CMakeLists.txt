include(CTest)

add_executable(test_material_composition test_material_composition.F90)
target_link_libraries(test_material_composition truchas-extra)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ./test_material_composition ${MPIEXEC_POSTFLAGS})
add_test(material_composition-1p ${test_launch})
set_tests_properties(material_composition-1p PROPERTIES PROCESSORS 1)
set_property(TEST material_composition-1p PROPERTY LABELS UNIT 2d-truchas)

set(test_launch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ./test_material_composition ${MPIEXEC_POSTFLAGS})
add_test(material_composition-4p ${test_launch})
set_tests_properties(material_composition-4p PROPERTIES PROCESSORS 4)
set_property(TEST material_composition-4p PROPERTY LABELS UNIT 2d-truchas)
