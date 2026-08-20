# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Truchas_MODULE_DIR})
include_directories(${PGSLib_MODULE_DIR})
include_directories(${TruchasExe_SOURCE_DIR}/utilities)

set(VOLUME_TRACKING_2D_SOURCE_FILES
  geom_axisymmetric.F90
  volume_tracker_2d_class.F90
  simple_volume_tracker_type.F90
  gradient_2d_cc_function.F90
  cell_geom_2d_vof_type.F90
  plane_2d_type.F90
  locate_plane_os_2d_function.F90
  truncation_volume_2d_type.F90
  geometric_volume_tracker_type.F90
)

add_library(volume_tracking_2d ${VOLUME_TRACKING_2D_SOURCE_FILES})
set_target_properties(volume_tracking_2d PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
target_link_libraries(volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_include_directories(volume_tracking_2d INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_library(test_vof_2d OBJECT read_inputfile.F90 vof_2d_test_driver.F90
            gaussian_quadrature_vofinit.F90)
target_link_libraries(test_vof_2d volume_tracking_2d)

add_executable(vof_advection vof_advection.F90 $<TARGET_OBJECTS:test_vof_2d>)
add_executable(vof_vortex vof_vortex.F90 $<TARGET_OBJECTS:test_vof_2d>)
add_executable(vof_axisymmetric vof_axisymmetric.F90 $<TARGET_OBJECTS:test_vof_2d>)
add_executable(vof_axialflow vof_axialflow.F90 $<TARGET_OBJECTS:test_vof_2d>)
add_executable(trunc_vol truncvolume_2dtest_driver.F90 $<TARGET_OBJECTS:test_vof_2d>)

target_link_libraries(vof_advection volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_link_libraries(vof_vortex volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_link_libraries(vof_axisymmetric volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_link_libraries(vof_axialflow volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})
target_link_libraries(trunc_vol volume_tracking_2d truchas ${TRUCHAS_LINK_LIBS})

configure_file(input_advection.txt input_advection.txt COPYONLY)
configure_file(input_vortex.txt input_vortex.txt COPYONLY)
configure_file(input_axisymmetric.txt input_axisymmetric.txt COPYONLY)
configure_file(input_axisymmetric_mixedmesh.txt input_axisymmetric_mixedmesh.txt COPYONLY)
configure_file(input_axialflow.txt input_axialflow.txt COPYONLY)
configure_file(circleadv_vof.txt circleadv_vof.txt COPYONLY)
configure_file(circlevort_vof.txt circlevort_vof.txt COPYONLY)
configure_file(circleaxisym_vof.txt circleaxisym_vof.txt COPYONLY)
configure_file(circleaxisym_vof_mixedmesh.txt circleaxisym_vof_mixedmesh.txt COPYONLY)
configure_file(circleaxialfl_vof.txt circleaxialfl_vof.txt COPYONLY)

add_test(truncated_volumes
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./trunc_vol ${MPIEXEC_POSTFLAGS})
set_tests_properties(truncated_volumes PROPERTIES PROCESSORS 1 LABELS 2d-truchas)

add_test(circle_advection
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./vof_advection ${MPIEXEC_POSTFLAGS})
set_tests_properties(circle_advection PROPERTIES PROCESSORS 4 LABELS 2d-truchas)

add_test(circle_vortex
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./vof_vortex ${MPIEXEC_POSTFLAGS})
set_tests_properties(circle_vortex PROPERTIES PROCESSORS 4 LABELS 2d-truchas)

add_test(circle_axisymm
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
    ./vof_axisymmetric input_axisymmetric.txt circleaxisym_vof.txt ${MPIEXEC_POSTFLAGS})
set_tests_properties(circle_axisymm PROPERTIES PROCESSORS 4 LABELS 2d-truchas)

add_test(circle_axisymm_mixedmesh
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
    ./vof_axisymmetric input_axisymmetric_mixedmesh.txt circleaxisym_vof_mixedmesh.txt ${MPIEXEC_POSTFLAGS})
set_tests_properties(circle_axisymm_mixedmesh PROPERTIES PROCESSORS 4 LABELS 2d-truchas)

add_test(circle_axialfl
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
    ./vof_axialflow input_axialflow.txt circleaxialfl_vof.txt ${MPIEXEC_POSTFLAGS})
set_tests_properties(circle_axialfl PROPERTIES PROCESSORS 4 LABELS 2d-truchas)
