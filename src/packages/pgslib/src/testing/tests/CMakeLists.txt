project(PGSLibTestMaster)

include(FortranPreprocessFiles)

include_directories(${PGSLib_MODULE_DIR})

set(MOD_F_FILES
        test-redux.F90
        test-bcast.F90
        test-dist.F90
        test-collate.F90
        error_tests.F90)
fortran_preprocess_files(MOD_SOURCE_FILES
                         FILES ${MOD_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestMaster_SOURCE_DIR}
			 PROCESS_TARGET ProcessTestMasterFiles)

fortran_preprocess_files(MAIN_SOURCE_FILE
                         FILES master-test.F90
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestMaster_SOURCE_DIR}
			 PROCESS_TARGET ProcessTestMasterMain)

# Add module library
add_library(testmaster ${MOD_SOURCE_FILES})
target_link_libraries(testmaster pgslib_fe)
add_dependencies(testmaster ProcessTestMasterFiles)

# MPI
if (ENABLE_MPI)
  add_executable(master-test ${MAIN_SOURCE_FILE})
  target_link_libraries(master-test testmaster pgslib_fe pgslib_mpi)
  add_dependencies(master-test ProcessTestMasterMain)
  add_test(master-mpi 
           ${PGSLib_PARALLEL_SCRIPT} master-test)
endif()

# Serial
add_executable(master-test-serial ${MAIN_SOURCE_FILE})
target_link_libraries(master-test-serial testmaster pgslib_fe pgslib_serial)
add_dependencies(master-test-serial ProcessTestMasterMain)
add_test(master-serial
         ${PGSLib_SERIAL_SCRIPT} master-test-serial)
