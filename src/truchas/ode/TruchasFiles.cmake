# Truchas files in directory
#   ode

# List of files to  process
set(ODE_FILES)

# List of files to add to the Truchas library
set(ODE_SOURCE_FILES)

# Process target name
set(ODE_TARGET_NAME ProcessTruchasOdeFiles)


set(ODE_FILES
       ode/bdf2/bdf2_controller.F90
       ode/bdf2/bdf2_dae.F90
       ode/bdf2/bdf2_integrator.F90
       ode/bdf2/bdf2_kinds.F90
       ode/bdf2/bdf2_profiling.F90
       ode/bdf2/solution_history.F90)

set(ODE_FPP_FLAGS 
        -I${TruchasExe_SOURCE_DIR}/ode
	-I${TruchasExe_SOURCE_DIR}/utilities
	${Truchas_FPP_FLAGS})

# Process files
fortran_preprocess_files(ODE_SOURCE_FILES
                         FILES ${ODE_FILES}
			 FPP_EXECUTABLE ${Truchas_PREPROCESSOR}
			 FPP_FLAGS ${ODE_FPP_FLAGS}
			 PROCESS_TARGET ${ODE_TARGET_NAME})
set_source_files_properties(${ODE_SOURCE_FILES} PROPERTIES
                            COMPILE_FLAGS -I${PGSLib_MODULE_DIR})

list(APPEND Truchas_LIBRARY_SOURCE_FILES ${ODE_SOURCE_FILES})		       
list(APPEND Truchas_PROCESS_TARGETS ${ODE_TARGET_NAME})


