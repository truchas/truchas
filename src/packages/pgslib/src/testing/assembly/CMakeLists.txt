project(PGSLibAssembly)

include(FortranPreprocessFiles)

set(MOD_F_FILES
        test_scan.F90
	test_sort.F90)
fortran_preprocess_files(MOD_SOURCE_FILES
                         FILES ${MOD_F_FILES}
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibAssembly_SOURCE_DIR}
			 PROCESS_TARGET ProcessAssemblyFiles)

fortran_preprocess_files(MAIN_SOURCE_FILE
                         FILES test-assembly.F90
			 GREP_FILTERS !!CPP!!
			 FPP_EXECUTABLE ${PGSLib_PREPROCESSOR}
			 FPP_FLAGS -I${PGSLibAssembly_SOURCE_DIR}
			 PROCESS_TARGET ProcessAssemblyMain)

# Add module library
add_library(testassem ${MOD_SOURCE_FILES})
add_dependencies(testassem pgslib_fe ProcessAssemblyFiles)

# MPI 	 
if (ENABLE_MPI)
  add_executable(test-assembly ${MAIN_SOURCE_FILE})
  target_link_libraries(test-assembly testassem pgslib_fe pgslib_mpi)
  add_dependencies(test-assembly ProcessAssemblyMain)
  add_test(assembly-mpi 
           ${PGSLib_PARALLEL_SCRIPT} test-assembly)
endif()	 

# Serial
add_executable(test-assembly-serial ${MAIN_SOURCE_FILE})
target_link_libraries(test-assembly-serial testassem pgslib_fe pgslib_serial)
add_dependencies(test-assembly-serial ProcessAssemblyMain)
add_test(assembly-serial 
         ${PGSLib_SERIAL_SCRIPT} test-assembly-serial)




	 
