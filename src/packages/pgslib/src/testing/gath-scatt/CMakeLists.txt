project(PGSLibTestGS)

include(FortranPreprocessFiles)

set(MOD_F_FILES
      gather_test.F90
      scatter_test.F90
      test_globals_module.F90
      test_scatter_module.F90
      test_gather_module.F90
      test_gs_setup_module.F90
      test_gs_parallel_module.F90)
fortran_preprocess_files(MOD_SOURCE_FILES
                         FILES ${MOD_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestGS_SOURCE_DIR}
			 PROCESS_TARGET ProcessGSTestFiles)

# Add a module library       
add_library(testgs ${MOD_SOURCE_FILES})
target_link_libraries(testgs pgslib_fe)

# MAIN (program) source file
set(MAIN_F_FILE test-gs.F90)
fortran_preprocess_files(MAIN_SOURCE_FILE
                         FILES ${MAIN_F_FILE}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestGS_SOURCE_DIR}
			 PROCESS_TARGET ProcessGSTestMain)

# MPI
if (ENABLE_MPI)
  add_executable(test-gs ${MAIN_SOURCE_FILE})
  target_link_libraries(test-gs testgs pgslib_fe pgslib_mpi)
  add_dependencies(test-gs ProcessGSTestFiles ProcessGSTestMain)
  add_test(gath-scatt-mpi 
           ${PGSLib_PARALLEL_SCRIPT} test-gs)
endif()

# Serial
add_executable(test-gs-serial ${MAIN_SOURCE_FILE})
target_link_libraries(test-gs-serial testgs pgslib_fe pgslib_serial)
add_dependencies(test-gs-serial ProcessGSTestFiles ProcessGSTestMain)
add_test(gath-scatt-serial 
         ${PGSLib_SERIAL_SCRIPT} test-gs-serial)




	 
