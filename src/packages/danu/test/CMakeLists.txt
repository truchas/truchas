# ############################################################################ #
#
#
#   DANU CMake 
#       Test Suite
#
# ############################################################################ #

include(DanuGlobalMacros)
include(PrintVariable)

if(ENABLE_Tests)

  project(DanuTestSuite)
  message(STATUS "Configure test suite")

  # Define the include directories
  include_directories(${HDF5_INCLUDE_DIRS})
  include_directories(${Check_INCLUDE_DIRS})

  include_directories("${Danu_SOURCE_DIR}/src")
  include_directories("${Danu_BINARY_DIR}/src/mod_files")

  # Test unit utilities (C interface utest)
  add_library(utest unit_test_utils.c random_generators.c)
  if ( HDF5_HL_FOUND )
    set(hdf5_libs ${HDF5_HL_LIBRARY} ${HDF5_C_LIBRARY} ${HDF5_LINK_LIBRARIES})
  else()
    set(hdf5_libs ${HDF5_C_LIBRARY} ${HDF5_LINK_LIBRARIES})
  endif()  
  target_link_libraries(utest danu ${hdf5_libs})
  if ( ${Check_BUILD_TARGET} )
    add_dependencies(utest ${Check_BUILD_TARGET})
  endif()

  # Test unit utilities (Fortran interface utest_f)
  if(ENABLE_Fortran)

    # Location of the *.mod files
    include_directories("${Danu_BINARY_DIR}/src")
    add_library(utest_f funit_utils.f90)
    target_link_Libraries(utest_f utest)
    add_library(fruit fruit.f90 fruit_util.f90)

  endif(ENABLE_Fortran)

  # Unit test utils
  set(UTEST_C_TESTS
              check_random.c)
  foreach(test ${UTEST_C_TESTS})
    #print_variable(test)
    add_check_test(${test} utest ${Check_LIBRARIES})
  endforeach()    
     
  # Danu C interface tests
  set(DANU_C_TESTS
           check_memory.c 
           check_valid.c
           check_file.c
           check_attribute.c 
           check_group.c
	   check_datasets.c
           check_output.c 
           check_mesh.c
           check_simulation.c
           check_non-series.c
           check_series.c
	   check_probes.c)

  # Need to remove the attribute test if HDF5_HL_LIBRARY is not found
  if ( NOT HDF5_HL_FOUND )
    list(REMOVE_ITEM DANU_C_TESTS check_attribute.c)
  endif() 

  foreach(test ${DANU_C_TESTS})
    #print_variable(test)
    add_check_test(${test} danu utest ${Check_LIBRARIES})
  endforeach() 

  if (ENABLE_Fortran)
    add_check_test(check_strings.c danu danu_f utest utest_f ${Check_LIBRARIES})
  endif()  

  # Fortran interface tests
  if(ENABLE_Fortran)
    set(FRUIT_TESTS
              funit_strings.f90
              funit_hidt.f90
              funit_file.f90
              funit_attribute.f90
              funit_group.f90
              funit_output.f90
              funit_mesh.f90
              funit_simulation.f90
              funit_non-series.f90
              funit_series-group.f90
              funit_series-data.f90
              funit_probes.f90)
    foreach(test ${FRUIT_TESTS})
      add_fruit_test(${test} danu_f utest_f fruit)
    endforeach()

  endif(ENABLE_Fortran)  



endif(ENABLE_Tests)  
