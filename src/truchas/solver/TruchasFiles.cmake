# Truchas files in directory
#   solver

# List of files to  process
set(SOLVER_FILES)

# List of files to add to the Truchas library
set(SOLVER_SOURCE_FILES)

# Process target name
set(SOLVER_TARGET_NAME ProcessTruchasSolverFiles)


set(SOLVER_FILES
          solver/hypre_c_binding.F90
          solver/fhypre.F90
          solver/hypre_pcg_type.F90
          solver/linear_solution.F90
          solver/nonlinear_solution.F90
          solver/ortho_matvec.F90
          solver/preconditioners.F90
          solver/ridders_type.F90)

set(SOLVER_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/utilities
        -I${TruchasExe_SOURCE_DIR}/solver
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(SOLVER_SOURCE_FILES
                         FILES ${SOLVER_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${SOLVER_FPP_FLAGS}
			 PROCESS_TARGET ${SOLVER_TARGET_NAME})
set(fc_flags)
if(ENABLE_PGSLib)
  list(APPEND fc_flags -I${PGSLib_MODULE_DIR})
endif()
if(ENABLE_UbikSolve)
  list(APPEND fc_flags -I${UbikSolve_MODULE_DIR})
endif()
build_whitespace_string(SOLVER_COMPILE_FLAGS ${fc_flags})
set(SOLVER_COMPILE_FLAGS "-I${PGSLib_MODULE_DIR} -I${UbikSolve_MODULE_DIR}")		       
set_source_files_properties(${SOLVER_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS ${SOLVER_COMPILE_FLAGS})

# HYPRE C File, needs compile flags                           
list(APPEND SOLVER_SOURCE_FILES solver/hypre_ext.c)

# Build compile flags string from list HYPRE_INCLUDE_DIRS and
# MPI_C_COMPILE_FLAGS
set(hypre_compile_flags)
set(hypre_flag_list)
foreach ( dir ${HYPRE_INCLUDE_DIRS} )
  list(APPEND hypre_flag_list "-I${dir}")
endforeach()
if(ENABLE_MPI AND DEFINED MPI_C_COMPILE_FLAGS)
  list(APPEND hypre_flag_list ${MPI_C_COMPILE_FLAGS})
endif()
build_whitespace_string(hypre_compile_flags ${hypre_flag_list})

set_source_files_properties(solver/hypre_ext.c PROPERTIES
                            COMPILE_FLAGS ${hypre_compile_flags})



list(APPEND Truchas_LIBRARY_SOURCE_FILES ${SOLVER_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${SOLVER_TARGET_NAME})


