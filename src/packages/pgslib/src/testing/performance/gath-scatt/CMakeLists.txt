project(PGSLibTestGSPerf)

include(FortranPreprocessFiles)

set(MOD_F_FILES
        perf_types.F90
        gather_perf_measure.F90)
fortran_preprocess_files(MOD_SOURCE_FILES
                         FILES ${MOD_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestGSPerf_SOURCE_DIR}
			 PROCESS_TARGET ProcessTestGSPerfFiles)

fortran_preprocess_files(MAIN_SOURCE_FILE
                         FILES gs-perf.F90
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibTestGSPerf_SOURCE_DIR}
			 PROCESS_TARGET ProcessTestGSPerfMain)

# Add a test library		       
add_library(testgsperf ${MOD_SOURCE_FILES})
target_link_libraries(testgsperf pgslib_fe)
add_dependencies(testgsperf ProcessTestGSPerfFiles)

# MPI
if (ENABLE_MPI)
  add_executable(gs-perf ${MAIN_SOURCE_FILE})
  target_link_libraries(gs-perf testgsperf pgslib_fe pgslib_mpi)
  add_dependencies(gs-perf ProcessTestGSPerfMain)
  set(MPI_ARGS ${MPIEXEC_NUMPROC_FLAG} 4 ${PGSLibTestGSPerf_BINARY_DIR}/gs-perf)
  add_test(gs-perf-mpi 
           ${MPIEXEC} ${MPI_ARGS})
endif()

# Serial
add_executable(gs-perf-serial ${MAIN_SOURCE_FILE})
target_link_libraries(gs-perf-serial testgsperf pgslib_fe pgslib_serial)
add_dependencies(gs-perf-serial ProcessTestGSPerfMain)
add_test(gs-perf-serial gs-perf-serial)
	 
